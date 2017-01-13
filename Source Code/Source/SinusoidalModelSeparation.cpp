#include "SinusoidalModelSeparation.h"
#include "Cluster.h"
#include "SinusoidalTrajectory.h"
#include "SinusoidalTrajectoryPoint.h"
#include "Transform.h"
#include <complex>
#include <limits>
#include <iostream>

//Helper method for identifying peak positions
//For a quadratic through points (-1,y_1),(0,y0),(1,y1), returns x-coordinate of peak/trough
inline double quadraticFitTurningPoint(double y_1, double y0, double y1)
{
	/*double ab = y1-y0;
	double a_b = y_1-y0;
	double a = (ab + a_b)/2;
	double b = (ab - a_b)/2;
	double grad = 2*a;
	double intercept = b;
	double root = -intercept/grad;*/
	return -0.5*(y1-y_1)/(y1+y_1-2*y0);
}

//Helper method for extracting sinusoids
//Given a vector representing samples from a continuous signal, gives approximate positions of peaks
vector<double>* peakPositions(vector<double>* values)
{
	vector<double>* retVal = new vector<double>();
	double prev = 0.;
	double current = (*values)[0];
	double next = (*values)[1];
	for(long i=2; i<values->size(); i++)
	{
		prev = current;
		current = next;
		next = (*values)[i];
		if(current > prev && current > next)
		{
			retVal->push_back(i-1+quadraticFitTurningPoint(prev, current, next));
		}
	}

	return retVal;
}

double approxPowerAroundPosition(vector<double>* values, double position)
{
	long lPos = (long) position;
	double part = position-lPos;
	return (1-part)*(*values)[lPos]*(*values)[lPos] + part*(*values)[lPos+1]*(*values)[lPos+1];
}

double approxPowerAroundPosition(vector<complex<double>>* values, double position)
{
	long lPos = (long) position;
	double part = position-lPos;
	return (1-part)*norm((*values)[lPos]) + part*norm((*values)[lPos+1]);
}

double stereoScaleAroundPosition(vector<complex<double>>* left, vector<complex<double>>* right, double position)
{
	double lPower = approxPowerAroundPosition(left, position);
	double rPower = approxPowerAroundPosition(right, position);
	return lPower / (lPower + rPower);
}

vector<vector<vector<double>>>* SinusoidalModelSeparation::separate(vector<double>* lSampleVector, vector<double>* rSampleVector, int numOfSources)
{
	//Obtain spectra
	vector<vector<complex<double>>>* lSpectrum = Transform::stft(lSampleVector, windowSize, hopSize, nfft);
	vector<vector<complex<double>>>* rSpectrum = Transform::stft(rSampleVector, windowSize, hopSize, nfft);

	//Get sizes
	long tCount = lSpectrum->size();
	long fCount = (*lSpectrum)[0].size();

	//Get overall power spectrum
	vector<vector<double>> absSpec = vector<vector<double>>();
	for(long t=0; t<tCount; t++)
	{
		absSpec.push_back(vector<double>());
		for(long f=0; f<fCount; f++)
		{
			absSpec[t].push_back(norm((*lSpectrum)[t][f]) + norm((*rSpectrum)[t][f]));
		}
	}
	
	cout << "Power spectrum obtained" << endl;

	//Estimate background noise level
	//NOTE::Could cluster power values into two clusters and take midpoint of centres as threshold
	//--Find row of power spectrum with lowest maximum
	double minRowMax = numeric_limits<double>::max();
	for(long f=0; f<fCount; f++)
	{
		double rowMax = numeric_limits<double>::min();
		for(long t=0; t<tCount; t++)
		{
			if(rowMax < absSpec[t][f])
			{
				rowMax = absSpec[t][f];
			}
		}
		if(minRowMax > rowMax)
		{
			minRowMax = rowMax;
		}
	}

	//--Calculate threshold to eliminate noise
	const double threshold = minRowMax * 2000000.;

	//--Threshold power spectrum
	for(long t=0; t<tCount; t++)
	{
		for(long f=0; f<fCount; f++)
		{
			if(threshold > absSpec[t][f])
			{
				absSpec[t][f] = 0;
			}
		}
	}

	cout << "Power spectrum thresholded" << endl;

	//Extract sinuosoids
	vector<SinusoidalTrajectory> sins = vector<SinusoidalTrajectory>();
	for(long t=0; t<tCount; t++)
	{
		//Find local peak frequencies for current time frame
		vector<double>* peakPos = peakPositions(&(absSpec[t]));
		for(long i=0; i<peakPos->size(); i++)
		{
			//Create a trajectory point for each peak
			double p = (*peakPos)[i];
			SinusoidalTrajectoryPoint point(approxPowerAroundPosition(&(absSpec[t]),p),
				p,
				stereoScaleAroundPosition(&((*lSpectrum)[t]), &((*rSpectrum)[t]), p));
			
			//Extend an existing trajectory if possible
			bool newSin = true;
			for(long i=0; i<sins.size(); i++)
			{
				if(sins[i].endIndex == t-1 && sins[i].canAccept(point))
				{
					sins[i].addPoint(point);
					newSin = false;
					break;
				}
			}

			//If no trajectory could be extended, create a new one
			if(newSin)
			{
				sins.push_back(SinusoidalTrajectory(t,point));
			}
		}

		delete peakPos;
	}

	cout << "Sinusoids extracted" << endl;

	//Normalise sinusoids
	long numOfSins = sins.size();
	double minFreq = numeric_limits<double>::max();
	for(long i=0; i<numOfSins; i++)
	{
		sins[i].normalise();
		if(minFreq > sins[i].meanFreq)
		{
			minFreq = sins[i].meanFreq;
		}
	}

	cout << "Sinusoids normalised" << endl;

	//Calculate distance matrix
	vector<vector<double>> d = vector<vector<double>>();
	for(long i=0; i<numOfSins; i++)
	{
		vector<double> di = vector<double>();

		//Distance matrix should be symmetric, so copy previously-calculated cells
		for(long j=0; j<i; j++)
		{
			di.push_back(d[j][i]);
		}

		//Calculate the new distances
		for(long j=i; j<numOfSins; j++)
		{
			di.push_back(SinusoidalTrajectory::distance(sins[i], sins[j], minFreq));
		}
		d.push_back(di);
	}

	cout << "Distance matrix calculated" << endl;

	//Cluster sinusoids
	vector<int>* clusterTags = Cluster::kClusterLloyd(&d, numOfSources);

	cout << "Sinusoids clustered" << endl;

	//Reconstruct sources
	vector<vector<vector<double>>>* retVal = new vector<vector<vector<double>>>();
	for(int s=0; s<numOfSources; s++)
	{
		//Initialise spectra as zero-matrices
		vector<vector<complex<double>>> lSoundSpec = vector<vector<complex<double>>>();
		vector<vector<complex<double>>> rSoundSpec = vector<vector<complex<double>>>();
		for(long t=0; t<tCount; t++)
		{
			lSoundSpec.push_back(vector<complex<double>>(fCount));
			rSoundSpec.push_back(vector<complex<double>>(fCount));
		}

		//Add peaks from trajectories in the current cluster
		for(long i=0; i<numOfSins; i++)
		{
			if((*clusterTags)[i] != s)
			{
				continue;
			}
			for(long t=sins[i].startIndex; t<=sins[i].endIndex; t++)
			{
				long freq = (long) sins[i].data[t-sins[i].startIndex].trueFreq;
				lSoundSpec[t][freq] = (*lSpectrum)[t][freq];
				lSoundSpec[t][freq+1] = (*lSpectrum)[t][freq+1];
				rSoundSpec[t][freq] = (*rSpectrum)[t][freq];
				rSoundSpec[t][freq+1] = (*rSpectrum)[t][freq+1];
			}
		}

		//Calculate the samples from the spectra
		vector<double>* lSamples = Transform::istft(&lSoundSpec, hopSize);
		vector<double>* rSamples = Transform::istft(&rSoundSpec, hopSize);
		retVal->push_back(vector<vector<double>>());
		retVal->back().push_back(vector<double>(*lSamples));
		retVal->back().push_back(vector<double>(*rSamples));

		delete lSamples;
		delete rSamples;
	}

	//Clean up
	delete lSpectrum;
	delete rSpectrum;
	delete clusterTags;

	return retVal;
}
