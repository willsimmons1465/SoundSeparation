#include "NMFSeparation.h"
#include "Matrix.h"
#include "Transform.h"
#include <iostream>

vector<vector<vector<double>>>* NMFSeparation::separate(vector<double>* lSampleVector, vector<double>* rSampleVector, int numOfSources)
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
	vector<Matrix<double>>* nmfResults = Transform::nmf(&absSpec, numOfSources);
	Matrix<double> mix = nmfResults->data()[0];
	Matrix<double> source = nmfResults->data()[1];

	cout << "NMF computed" << endl;

	//Reconstruct each source
	vector<vector<vector<double>>>* retVal = new vector<vector<vector<double>>>();
	for(int s=0; s<numOfSources; s++)
	{
		//Obtain absolute spectrum for source
		Matrix<double> soundMix = Matrix<double>(tCount,1);
		Matrix<double> soundSource = Matrix<double>(1,fCount);
		for(long t=0; t<tCount; t++)
		{
			soundMix.data[t][0] = mix.data[t][s];
		}
		for(long f=0; f<fCount; f++)
		{
			soundSource.data[0][f] = source.data[s][f];
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