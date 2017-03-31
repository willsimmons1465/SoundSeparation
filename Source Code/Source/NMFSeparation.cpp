#include "NMFSeparation.h"
#include "Cluster.h"
//#include "Matrix.h"
//#include <Eigen/Core>
#include "Transform.h"
#include <iostream>

vector<vector<vector<double>>>* NMFSeparation::separate(vector<double>* lSampleVector, vector<double>* rSampleVector, int numOfSources, int factors)
{
	//Obtain spectra
	Eigen::MatrixXcd* lSpectrum = Transform::stft(lSampleVector, windowSize, hopSize, nfft);
	Eigen::MatrixXcd* rSpectrum = Transform::stft(rSampleVector, windowSize, hopSize, nfft);

	cout << "Spectra obtained" << endl;

	//Get sizes
	long tCount = lSpectrum->rows();
	long fCount = lSpectrum->cols();

	//Get overall power spectrum
	Eigen::MatrixXd absSpec = Eigen::MatrixXd(tCount,fCount);
	for(long t=0; t<tCount; t++)
	{
		for(long f=0; f<fCount; f++)
		{
			absSpec(t, f) = norm((*lSpectrum)(t, f)) + norm((*rSpectrum)(t, f));
		}
	}

	cout << "Power spectrum obtained" << endl;

	//Perform NMF
	vector<Eigen::MatrixXd>* nmfResults = Transform::nmf(&absSpec, factors);
	Eigen::MatrixXd mix = nmfResults->data()[0];
	Eigen::MatrixXd source = nmfResults->data()[1];

	cout << "NMF computed" << endl;

	//Construct feature vectors
	Eigen::MatrixXd features = Eigen::MatrixXd(factors, factors*2);
	for(int i=0; i<factors; i++)
	{
		for(int j=0; j<i; j++)
		{
			features(i, j) = features(j, i);
			features(i, factors+j) = features(j, factors+i);
		}
		for(int j=i; j<factors; j++)
		{
			double mixDistance = 0.;
			for(long t=0; t<tCount; t++)
			{
				mixDistance += pow(mix(t, i) - mix(t, j),2);
			}
			features(i, j) = mixDistance;

			double sourceDistance = 0.;
			for(long f=0; f<fCount; f++)
			{
				sourceDistance += pow(source(i, f) - source(j, f),2);
			}
			features(i, factors+j) = sourceDistance;
		}
	}

	cout << "Feature vectors created" << endl;

	//Cluster factors
	vector<int>* clusterTags = Cluster::kClusterLloyd(&features, numOfSources);

	cout << "Factors clustered" << endl;

	//Reconstruct each source
	vector<vector<vector<double>>>* retVal = new vector<vector<vector<double>>>();
	for(int s=0; s<numOfSources; s++)
	{
		//Obtain absolute spectrum for source
		Eigen::MatrixXd soundMix = Eigen::MatrixXd(tCount,0);
		Eigen::MatrixXd soundSource = Eigen::MatrixXd(0,fCount);
		for(int a=0; a<factors; a++)
		{
			if(clusterTags->data()[a] != s)
			{
				continue;
			}
			
			soundMix.conservativeResize(Eigen::NoChange, soundMix.cols()+1);
			soundMix.col(soundMix.cols()-1) = mix.col(a);
			soundSource.conservativeResize(soundSource.rows()+1, Eigen::NoChange);
			soundSource.row(soundSource.rows()-1) = source.row(a);
		}
		Eigen::MatrixXd soundAbsSpec = soundMix * soundSource;

		//Consider stereo balance and phase of original source.
		Eigen::MatrixXcd lSoundSpec = Eigen::MatrixXcd(tCount, fCount);
		Eigen::MatrixXcd rSoundSpec = Eigen::MatrixXcd(tCount, fCount);
		for(long t=0; t<tCount; t++)
		{
			for(long f=0; f<fCount; f++)
			{
				double stereo = norm((*lSpectrum)(t, f)) / absSpec(t, f);
				lSoundSpec(t, f) = polar(sqrt(stereo * soundAbsSpec(t, f)), arg((*lSpectrum)(t, f)));
				rSoundSpec(t, f) = polar(sqrt((1-stereo) * soundAbsSpec(t, f)), arg((*rSpectrum)(t, f)));
			}
		}

		//Calculate the samples from the spectra
		Eigen::VectorXd* lVector = Transform::istft(&lSoundSpec, hopSize);
		Eigen::VectorXd* rVector = Transform::istft(&rSoundSpec, hopSize);
		vector<double> lSamples = vector<double>(lVector->data(), lVector->data() + lVector->size());
		vector<double> rSamples = vector<double>(rVector->data(), rVector->data() + rVector->size());
		retVal->push_back(vector<vector<double>>());
		retVal->back().push_back(vector<double>(lSamples));
		retVal->back().push_back(vector<double>(rSamples));

		delete lVector;
		delete rVector;
	}

	//Clean up
	delete lSpectrum;
	delete rSpectrum;
	delete nmfResults;

	return retVal;
}