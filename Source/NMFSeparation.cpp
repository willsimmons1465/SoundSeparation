#include "NMFSeparation.h"
#include "Cluster.h"
#include "Matrix.h"
#include "Transform.h"
#include <iostream>

vector<vector<vector<double>>>* NMFSeparation::separate(vector<double>* lSampleVector, vector<double>* rSampleVector, int numOfSources, int factors)
{
	//Obtain spectra
	Matrix<complex<double>>* lSpectrum = Transform::stft(lSampleVector, windowSize, hopSize, nfft);
	Matrix<complex<double>>* rSpectrum = Transform::stft(rSampleVector, windowSize, hopSize, nfft);

	cout << "Spectra obtained" << endl;

	//Get sizes
	long tCount = lSpectrum->data.size();
	long fCount = lSpectrum->data[0].size();

	//Get overall power spectrum
	Matrix<double> absSpec = Matrix<double>(tCount,fCount);
	for(long t=0; t<tCount; t++)
	{
		for(long f=0; f<fCount; f++)
		{
			absSpec.data[t][f] = norm(lSpectrum->data[t][f]) + norm(rSpectrum->data[t][f]);
		}
	}

	cout << "Power spectrum obtained" << endl;

	//Perform NMF
	vector<Matrix<double>>* nmfResults = Transform::nmf(&absSpec, factors);
	Matrix<double> mix = nmfResults->data()[0];
	Matrix<double> source = nmfResults->data()[1];

	cout << "NMF computed" << endl;

	//Construct feature vectors
	Matrix<double> features = Matrix<double>(factors, factors*2);
	for(int i=0; i<factors; i++)
	{
		for(int j=0; j<i; j++)
		{
			features.data[i][j] = features.data[j][i];
			features.data[i][factors+j] = features.data[j][factors+i];
		}
		for(int j=i; j<factors; j++)
		{
			double mixDistance = 0.;
			for(long t=0; t<tCount; t++)
			{
				mixDistance += pow(mix.data[t][i] - mix.data[t][j],2);
			}
			features.data[i][j] = mixDistance;

			double sourceDistance = 0.;
			for(long f=0; f<fCount; f++)
			{
				sourceDistance += pow(source.data[i][f] - source.data[j][f],2);
			}
			features.data[i][factors+j] = sourceDistance;
		}
	}

	cout << "Feature vectors created" << endl;

	//Cluster factors
	vector<int>* clusterTags = Cluster::kClusterLloyd(&(features.data), numOfSources);

	cout << "Factors clustered" << endl;

	//Reconstruct each source
	vector<vector<vector<double>>>* retVal = new vector<vector<vector<double>>>();
	for(int s=0; s<numOfSources; s++)
	{
		//Obtain absolute spectrum for source
		Matrix<double> soundMix = Matrix<double>(tCount,0);
		Matrix<double> soundSource = Matrix<double>(0,fCount);
		for(int a=0; a<factors; a++)
		{
			if(clusterTags->data()[a] != s)
			{
				continue;
			}
			
			for(long t=0; t<tCount; t++)
			{
				soundMix.data[t].push_back(mix.data[t][a]);
			}
			soundSource.data.push_back(source.data[a]);
		}
		Matrix<double>* soundAbsSpec = Matrix<double>::multiply(&soundMix, &soundSource);

		//Consider stereo balance and phase of original source.
		Matrix<complex<double>> lSoundSpec = Matrix<complex<double>>(tCount, fCount);
		Matrix<complex<double>> rSoundSpec = Matrix<complex<double>>(tCount, fCount);
		for(long t=0; t<tCount; t++)
		{
			for(long f=0; f<fCount; f++)
			{
				double stereo = norm(lSpectrum->data[t][f]) / absSpec.data[t][f];
				lSoundSpec.data[t][f] = polar(sqrt(stereo * soundAbsSpec->data[t][f]), arg(lSpectrum->data[t][f]));
				rSoundSpec.data[t][f] = polar(sqrt((1-stereo) * soundAbsSpec->data[t][f]), arg(rSpectrum->data[t][f]));
			}
		}

		//Calculate the samples from the spectra
		vector<double>* lSamples = Transform::istft(&lSoundSpec, hopSize);
		vector<double>* rSamples = Transform::istft(&rSoundSpec, hopSize);
		retVal->push_back(vector<vector<double>>());
		retVal->back().push_back(vector<double>(*lSamples));
		retVal->back().push_back(vector<double>(*rSamples));

		delete soundAbsSpec;
		delete lSamples;
		delete rSamples;
	}

	//Clean up
	delete lSpectrum;
	delete rSpectrum;
	delete nmfResults;

	return retVal;
}