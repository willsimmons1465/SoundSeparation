#include "SinusoidalModelSeparation.h"
#include "Cluster.h"
#include "Matrix.h"
#include "SinusoidalTrajectory.h"
#include "SinusoidalTrajectoryPoint.h"
#include "Transform.h"
#include "WavFileManager.h"
#include <complex>
#include <limits>
#include <iostream>

double clusterDistThreshold = 3.2;
double softClusterStiffness = 0.00022;

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
	Matrix<complex<double>>* lSpectrum = Transform::stft(lSampleVector, windowSize, hopSize, nfft);
	Matrix<complex<double>>* rSpectrum = Transform::stft(rSampleVector, windowSize, hopSize, nfft);

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
	
	//cout << "Power spectrum obtained" << endl;

	//Estimate background noise level
	//NOTE::Could cluster power values into two clusters and take midpoint of centres as threshold
	//--Find row of power spectrum with lowest maximum
	double minRowMax = numeric_limits<double>::max();
	for(long f=0; f<fCount; f++)
	{
		double rowMax = numeric_limits<double>::min();
		for(long t=0; t<tCount; t++)
		{
			if(rowMax < absSpec.data[t][f])
			{
				rowMax = absSpec.data[t][f];
			}
		}
		if(minRowMax > rowMax)
		{
			minRowMax = rowMax;
		}
	}

	//--Calculate threshold to eliminate noise
	const double threshold = minRowMax * 2000000.;
	//const double threshold = minRowMax * 100000.;
	//const double threshold = minRowMax * 100.;

	//--Threshold power spectrum
	for(long t=0; t<tCount; t++)
	{
		for(long f=0; f<fCount; f++)
		{
			if(threshold > absSpec.data[t][f])
			{
				absSpec.data[t][f] = 0;
			}
		}
	}

	//cout << "Power spectrum thresholded" << endl;

	//Extract sinuosoids
	vector<SinusoidalTrajectory> sins = vector<SinusoidalTrajectory>();
	for(long t=0; t<tCount; t++)
	{
		//Find local peak frequencies for current time frame
		vector<double>* peakPos = peakPositions(&(absSpec.data[t]));
		for(long i=0; i<peakPos->size(); i++)
		{
			//Create a trajectory point for each peak
			double p = (*peakPos)[i];
			SinusoidalTrajectoryPoint point(approxPowerAroundPosition(&(absSpec.data[t]),p),
				p,
				stereoScaleAroundPosition(&(lSpectrum->data[t]), &(rSpectrum->data[t]), p));
			
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

	//cout << "Sinusoids extracted" << endl;

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

	//cout << "Sinusoids normalised" << endl;

	//Calculate distance matrix
	/*Matrix<double> d = Matrix<double>(numOfSins, numOfSins);
	for(long i=0; i<numOfSins; i++)
	{
		//Distance matrix should be symmetric, so copy previously-calculated cells
		for(long j=0; j<i; j++)
		{
			d.data[i][j] = d.data[j][i];
		}

		//Calculate the new distances
		for(long j=i; j<numOfSins; j++)
		{
			d.data[i][j] = SinusoidalTrajectory::distance(sins[i], sins[j], minFreq);
		}
	}*/
	Matrix<double> d = Matrix<double>(numOfSins, numOfSins*6);
	for(long i=0; i<numOfSins; i++)
	{
		for(long j=i; j<numOfSins; j++)
		{
			d.data[i][6*j] = SinusoidalTrajectory::distanceHarm(sins[i], sins[j], minFreq);
			d.data[i][6*j+1] = SinusoidalTrajectory::distanceStereo(sins[i], sins[j]);
			d.data[i][6*j+2] = SinusoidalTrajectory::distanceOnset(sins[i], sins[j]);
			d.data[i][6*j+3] = SinusoidalTrajectory::distanceMiss(sins[i], sins[j]);
			d.data[i][6*j+4] = SinusoidalTrajectory::distanceAmp(sins[i], sins[j]);
			d.data[i][6*j+5] = SinusoidalTrajectory::distanceFreq(sins[i], sins[j]);

			if(j!=i)
			{
				d.data[j][6*i] = d.data[i][6*j];
				d.data[j][6*i+1] = d.data[i][6*j+1];
				d.data[j][6*i+2] = d.data[i][6*j+2];
				d.data[j][6*i+3] = d.data[i][6*j+3];
				d.data[j][6*i+4] = d.data[i][6*j+4];
				d.data[j][6*i+5] = d.data[i][6*j+5];
			}
		}
	}

	//cout << "Distance matrix calculated" << endl;

	//Cluster sinusoids
	//vector<int>* clusterTags = Cluster::kClusterLloyd(&(d.data), numOfSources);

	/*vector<int>* clusterTags = new vector<int>(numOfSins);
	for(long i=0; i<numOfSins; i++)
	{
		(*clusterTags)[i] = -1;
	}
	for(int s=0; s<numOfSources; s++)
	{
		long seed = -1;
		double maxAmp = 0.;
		for(long i=0; i<numOfSins; i++)
		{
			if((*clusterTags)[i] == -1 && maxAmp < sins[i].meanAmp)
			{
				seed = i;
				maxAmp = sins[i].meanAmp;
			}
		}
		if (seed == -1)
		{
			break;
		}
		for(long i=0; i<numOfSins; i++)
		{
			if(SinusoidalTrajectory::distanceHarm(sins[i], sins[seed], minFreq) < clusterDistThreshold)
			{
				(*clusterTags)[i] = s;
			}
		}
	}
	for(long i=0; i<numOfSins; i++)
	{
		if((*clusterTags)[i] != -1)
		{
			continue;
		}
		vector<double> maxHarmDists = vector<double>(numOfSources);
		for(long j=0; j<numOfSins; j++)
		{
			int jTag = (*clusterTags)[j];
			if(jTag != -1)
			{
				double dh = SinusoidalTrajectory::distanceHarm(sins[i], sins[j], minFreq);
				if(maxHarmDists[jTag] < dh)
				{
					maxHarmDists[jTag] = dh;
				}
			}
		}
		double minMaxHarmDist = numeric_limits<double>::max();
		for(int s=0; s<numOfSources; s++)
		{
			if(minMaxHarmDist > maxHarmDists[s])
			{
				minMaxHarmDist = maxHarmDists[s];
				(*clusterTags)[i] = s;
			}
		}
	}*/

	//Matrix<double>* clusterWeights = Cluster::softKCluster(&(d.data), numOfSources, 0.00001);

	Matrix<double>* clusterWeights = Cluster::softKCluster(&(d.data), numOfSources, softClusterStiffness);
	double max = 0;
	for (int i=0; i<numOfSins; i++)
		{
			for (int s=0; s<numOfSources; s++)
			{
				max += clusterWeights->data[i][s] * clusterWeights->data[i][s];
			}
		}
	for (int r=0; r<10; r++)
	{
		Matrix<double>* clusterWeightsTry = Cluster::softKCluster(&(d.data), numOfSources, softClusterStiffness);
		double sum = 0.;
		for (int i=0; i<numOfSins; i++)
		{
			for (int s=0; s<numOfSources; s++)
			{
				sum += clusterWeightsTry->data[i][s] * clusterWeightsTry->data[i][s];
			}
		}
		if (max < sum)
		{
			max = sum;
			delete clusterWeights;
			clusterWeights = clusterWeightsTry;
		}
		else
		{
			delete clusterWeightsTry;
		}
	}

	/*vector<Matrix<double>>* nmfResults = Transform::nmf(&d, numOfSources);

	for(long i=0; i<numOfSins; i++)
	{
		double sum = 0.;
		for(int s=0; s<numOfSources; s++)
		{
			sum += nmfResults->at(0).data[i][s];
		}
		if(sum == 0)
		{
			sum = 1e-307;
		}
		for(int s=0; s<numOfSources; s++)
		{
			nmfResults->at(0).data[i][s] /= sum;
		}
	}*/

	//cout << "Sinusoids clustered" << endl;

	//Reconstruct sources
	/*vector<vector<vector<double>>>* retVal = new vector<vector<vector<double>>>();
	for(int s=0; s<numOfSources; s++)
	{
		//Initialise spectra as zero-matrices
		Matrix<complex<double>> lSoundSpec = Matrix<complex<double>>(tCount, fCount);
		Matrix<complex<double>> rSoundSpec = Matrix<complex<double>>(tCount, fCount);

		//Add peaks from trajectories in the current cluster
		for(long i=0; i<numOfSins; i++)
		{
			/*if((*clusterTags)[i] != s)
			{
				continue;
			}*/
			/*for(long t=sins[i].startIndex; t<=sins[i].endIndex; t++)
			{
				long freq = (long) sins[i].data[t-sins[i].startIndex].trueFreq;
				/*lSoundSpec.data[t][freq] = lSpectrum->data[t][freq];
				lSoundSpec.data[t][freq+1] = lSpectrum->data[t][freq+1];
				rSoundSpec.data[t][freq] = rSpectrum->data[t][freq];
				rSoundSpec.data[t][freq+1] = rSpectrum->data[t][freq+1];*/
				/*lSoundSpec.data[t][freq] = lSpectrum->data[t][freq] * clusterWeights->data[i][s];
				lSoundSpec.data[t][freq+1] = lSpectrum->data[t][freq+1] * clusterWeights->data[i][s];
				rSoundSpec.data[t][freq] = rSpectrum->data[t][freq] * clusterWeights->data[i][s];
				rSoundSpec.data[t][freq+1] = rSpectrum->data[t][freq+1] * clusterWeights->data[i][s];*/
				/*lSoundSpec.data[t][freq] = lSpectrum->data[t][freq] * nmfResults->at(0).data[i][s];
				lSoundSpec.data[t][freq+1] = lSpectrum->data[t][freq+1] * nmfResults->at(0).data[i][s];
				rSoundSpec.data[t][freq] = rSpectrum->data[t][freq] * nmfResults->at(0).data[i][s];
				rSoundSpec.data[t][freq+1] = rSpectrum->data[t][freq+1] * nmfResults->at(0).data[i][s];
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
	}*/

	vector<vector<vector<double>>>* retVal = new vector<vector<vector<double>>>();
	vector<Matrix<complex<double>>> lSoundSpecs = vector<Matrix<complex<double>>>();
	vector<Matrix<complex<double>>> rSoundSpecs = vector<Matrix<complex<double>>>();
	for(int s=0; s<numOfSources; s++)
	{
		//Initialise spectra as zero-matrices
		lSoundSpecs.push_back(Matrix<complex<double>>(tCount, fCount));
		rSoundSpecs.push_back(Matrix<complex<double>>(tCount, fCount));

		//Add peaks from trajectories in the current cluster
		for(long i=0; i<numOfSins; i++)
		{
			for(long t=sins[i].startIndex; t<=sins[i].endIndex; t++)
			{
				long freq = (long) sins[i].data[t-sins[i].startIndex].trueFreq;
				lSoundSpecs[s].data[t][freq] = lSpectrum->data[t][freq] * clusterWeights->data[i][s];
				lSoundSpecs[s].data[t][freq+1] = lSpectrum->data[t][freq+1] * clusterWeights->data[i][s];
				rSoundSpecs[s].data[t][freq] = rSpectrum->data[t][freq] * clusterWeights->data[i][s];
				rSoundSpecs[s].data[t][freq+1] = rSpectrum->data[t][freq+1] * clusterWeights->data[i][s];
				/*lSoundSpecs[s].data[t][freq] = lSpectrum->data[t][freq] * nmfResults->at(0).data[i][s];
				lSoundSpecs[s].data[t][freq+1] = lSpectrum->data[t][freq+1] * nmfResults->at(0).data[i][s];
				rSoundSpecs[s].data[t][freq] = rSpectrum->data[t][freq] * nmfResults->at(0).data[i][s];
				rSoundSpecs[s].data[t][freq+1] = rSpectrum->data[t][freq+1] * nmfResults->at(0).data[i][s];*/
			}
		}

		for(long t=0; t<tCount; t++)
		{
			for(long f=0; f<fCount; f++)
			{
				lSpectrum->data[t][f] -= lSoundSpecs[s].data[t][f];
				rSpectrum->data[t][f] -= rSoundSpecs[s].data[t][f];
			}
		}
	}

	for(int s=0; s<numOfSources; s++)
	{
		for(long t=0; t<tCount; t++)
		{
			for(long f=0; f<fCount; f++)
			{
				lSoundSpecs[s].data[t][f] += lSpectrum->data[t][f] / (double) numOfSources;
				rSoundSpecs[s].data[t][f] += rSpectrum->data[t][f] / (double) numOfSources;
			}
		}

		//Calculate the samples from the spectra
		vector<double>* lSamples = Transform::istft(&(lSoundSpecs[s]), hopSize);
		vector<double>* rSamples = Transform::istft(&(rSoundSpecs[s]), hopSize);
		retVal->push_back(vector<vector<double>>());
		retVal->back().push_back(vector<double>(*lSamples));
		retVal->back().push_back(vector<double>(*rSamples));

		delete lSamples;
		delete rSamples;
	}

	//Clean up
	delete lSpectrum;
	delete rSpectrum;
	//delete clusterTags;
	delete clusterWeights;
	//delete nmfResults;

	return retVal;
}

double euclideanSquaredDistance(vector<double>* a, vector<double>* b)
{
	vector<double>* smaller = (a->size() < b->size())?a:b;
	vector<double>* longer = (a->size() < b->size())?b:a;
	double retVal = 0.;
	for(long i=0; i<smaller->size(); i++)
	{
		retVal += pow(smaller->data()[i] - longer->data()[i], 2);
	}
	for(long i=smaller->size(); i<longer->size(); i++)
	{
		retVal += pow(longer->data()[i], 2);
	}

	return retVal;
}

double proportionalDifference(vector<double>* al, vector<double>* ar, vector<double>* bl, vector<double>* br)
{
	vector<double>* smallerl = (al->size() < bl->size())?al:bl;
	vector<double>* smallerr = (al->size() < bl->size())?ar:br;
	vector<double>* longerl = (al->size() < bl->size())?bl:al;
	vector<double>* longerr = (al->size() < bl->size())?br:ar;
	double rms1 = 0.;
	for(long i=0; i<smallerl->size(); i++)
	{
		rms1 += pow(smallerl->data()[i], 2);
		rms1 += pow(smallerr->data()[i], 2);
	}
	rms1 = sqrt(rms1 / smallerl->size());
	double rms2 = 0.;
	for(long i=0; i<longerl->size(); i++)
	{
		rms2 += pow(longerl->data()[i], 2);
		rms2 += pow(longerr->data()[i], 2);
	}
	rms2 = sqrt(rms2 / longerl->size());
	double diff = 0.;
	for(long i=0; i<smallerl->size(); i++)
	{
		diff += pow(smallerl->data()[i]/rms1 - longerl->data()[i]/rms2, 2);
		diff += pow(smallerr->data()[i]/rms1 - longerr->data()[i]/rms2, 2);
	}
	for(long i=smallerl->size(); i<longerl->size(); i++)
	{
		diff += pow(longerl->data()[i]/rms2, 2);
		diff += pow(longerr->data()[i]/rms2, 2);
	}
	
	return diff;
}

double stereoCosineSimilarity(vector<double>* al, vector<double>* ar, vector<double>* bl, vector<double>* br)
{
	vector<double>* smallerl = (al->size() < bl->size())?al:bl;
	vector<double>* smallerr = (al->size() < bl->size())?ar:br;
	vector<double>* longerl = (al->size() < bl->size())?bl:al;
	vector<double>* longerr = (al->size() < bl->size())?br:ar;
	double smallerLength = 0.;
	for(int i=0; i<smallerl->size(); i++)
	{
		smallerLength += pow(smallerl->data()[i], 2);
		smallerLength += pow(smallerr->data()[i], 2);
	}
	smallerLength = sqrt(smallerLength);
	double longerLength = 0.;
	for(int i=0; i<longerl->size(); i++)
	{
		longerLength += pow(longerl->data()[i], 2);
		longerLength += pow(longerr->data()[i], 2);
	}
	longerLength = sqrt(longerLength);
	double product = 0.;
	for(long i=0; i<smallerl->size(); i++)
	{
		product += smallerl->data()[i] * longerl->data()[i];
		product += smallerr->data()[i] * longerr->data()[i];
	}
	/*for(long i=smallerl->size(); i<longerl->size(); i++)
	{
		product += longerl->data()[i];
		product += longerr->data()[i];
	}*/
	return product / (smallerLength * longerLength);
}

double totalErrorFromTestSet(void)
{
	srand(19873495);

	//Get samples
	WavFileManager guitar("OptimiseTests/guitar_D3_very-long_forte_normal");
	WavFileManager trumpet("OptimiseTests/trumpet_Gs4_1_forte_normal");
	WavFileManager clarinet("OptimiseTests/clarinet_As3_1_forte_normal");
	WavFileManager saxophone("OptimiseTests/saxophone_E5_1_forte_normal");
	WavFileManager violin("OptimiseTests/violin_Fs6_1_forte_arco-normal");
	vector<WavFileManager*> sampleFiles = vector<WavFileManager*>();
	sampleFiles.push_back(&guitar);
	sampleFiles.push_back(&trumpet);
	sampleFiles.push_back(&clarinet);
	sampleFiles.push_back(&saxophone);
	sampleFiles.push_back(&violin);
	/*sampleFiles.push_back(WavFileManager("OptimiseTests/guitar_D3_very-long_forte_normal"));
	sampleFiles.push_back(WavFileManager("OptimiseTests/trumpet_Gs4_1_forte_normal"));
	sampleFiles.push_back(WavFileManager("OptimiseTests/clarinet_As3_1_forte_normal"));
	sampleFiles.push_back(WavFileManager("OptimiseTests/saxophone_E5_1_forte_normal"));
	sampleFiles.push_back(WavFileManager("OptimiseTests/violin_Fs6_1_forte_arco-normal"));*/

	//Find sum of distances between pairs of originals and separations
	double sumDistance = 0.;
	for(int i=0; i<sampleFiles.size(); i++)
	{
		for(int j=i+1; j<sampleFiles.size(); j++)
		{
			//Get samples from each sound
			vector<double> ilSamples;
			vector<double> irSamples;
			vector<double> jlSamples;
			vector<double> jrSamples;
			sampleFiles[i]->readSoundSample(ilSamples, irSamples);
			sampleFiles[j]->readSoundSample(jlSamples, jrSamples);

			//Mix sounds
			vector<double> mixLSamples;
			vector<double> mixRSamples;
			if(ilSamples.size() < jlSamples.size())
			{
				mixLSamples = vector<double>(jlSamples);
				mixRSamples = vector<double>(jrSamples);
				for(long k=0; k<ilSamples.size(); k++)
				{
					mixLSamples[k] += ilSamples[k];
					mixRSamples[k] += irSamples[k];
				}
			}
			else
			{
				mixLSamples = vector<double>(ilSamples);
				mixRSamples = vector<double>(irSamples);
				for(long k=0; k<jlSamples.size(); k++)
				{
					mixLSamples[k] += jlSamples[k];
					mixRSamples[k] += jrSamples[k];
				}
			}

			//Run separation
			vector<vector<vector<double>>>* separated = SinusoidalModelSeparation::separate(&mixLSamples, &mixRSamples, 2);

			//Calculate distances between results and originals
			double dist = euclideanSquaredDistance(&(separated->data()[0][0]), &ilSamples);
			dist += euclideanSquaredDistance(&(separated->data()[0][1]), &irSamples);
			dist += euclideanSquaredDistance(&(separated->data()[1][0]), &jlSamples);
			dist += euclideanSquaredDistance(&(separated->data()[1][1]), &jrSamples);
			double distSwitched = euclideanSquaredDistance(&(separated->data()[0][0]), &jlSamples);
			distSwitched += euclideanSquaredDistance(&(separated->data()[0][1]), &jrSamples);
			distSwitched += euclideanSquaredDistance(&(separated->data()[1][0]), &ilSamples);
			distSwitched += euclideanSquaredDistance(&(separated->data()[1][1]), &irSamples);
			sumDistance += (dist<distSwitched)?dist:distSwitched;

			cout << "*";

			delete separated;
		}
	}
	cout << endl;

	return sumDistance;
}

double totalErrorFromStereoTestSet(void)
{
	srand(19873495);

	//Get samples
	WavFileManager guitar("OptimiseTests/guitar_D3_very-long_forte_normal");
	WavFileManager trumpet("OptimiseTests/trumpet_Gs4_1_forte_normal");
	WavFileManager clarinet("OptimiseTests/clarinet_As3_1_forte_normal");
	WavFileManager saxophone("OptimiseTests/saxophone_E5_1_forte_normal");
	WavFileManager violin("OptimiseTests/violin_Fs6_1_forte_arco-normal");
	vector<WavFileManager*> sampleFiles = vector<WavFileManager*>();
	sampleFiles.push_back(&guitar);
	sampleFiles.push_back(&trumpet);
	sampleFiles.push_back(&clarinet);
	sampleFiles.push_back(&saxophone);
	sampleFiles.push_back(&violin);

	//Find sum of distances between pairs of originals and separations
	double sumDistance = 0.;
	for(int i=0; i<sampleFiles.size(); i++)
	{
		for(int j=i+1; j<sampleFiles.size(); j++)
		{
			//Get samples from each sound
			vector<double> ilSamples;
			vector<double> irSamples;
			vector<double> jlSamples;
			vector<double> jrSamples;
			sampleFiles[i]->readSoundSample(ilSamples, irSamples);
			sampleFiles[j]->readSoundSample(jlSamples, jrSamples);

			//Generate stereo positions
			double iStereo = ((i+j) % sampleFiles.size()) / (sampleFiles.size()-1);
			double jStereo = ((i-j) % sampleFiles.size()) / (sampleFiles.size()-1);

			//Mix sounds
			vector<double> mixLSamples;
			vector<double> mixRSamples;
			if(ilSamples.size() < jlSamples.size())
			{
				mixLSamples = vector<double>(jlSamples.size());
				mixRSamples = vector<double>(jrSamples.size());
			}
			else
			{
				mixLSamples = vector<double>(ilSamples.size());
				mixRSamples = vector<double>(irSamples.size());
			}
			for(long k=0; k<ilSamples.size(); k++)
			{
				ilSamples[k] = sqrt(1-iStereo)*(ilSamples[k] + irSamples[k])/2;
				irSamples[k] = sqrt(iStereo)*(ilSamples[k] + irSamples[k])/2;
				mixLSamples[k] += ilSamples[k];
				mixRSamples[k] += irSamples[k];
			}
			for(long k=0; k<jlSamples.size(); k++)
			{
				jlSamples[k] = sqrt(1-jStereo)*(jlSamples[k] + jrSamples[k])/2;
				jrSamples[k] = sqrt(jStereo)*(jlSamples[k] + jrSamples[k])/2;
				mixLSamples[k] += jlSamples[k];
				mixRSamples[k] += jrSamples[k];
			}

			//Run separation
			vector<vector<vector<double>>>* separated = SinusoidalModelSeparation::separate(&mixLSamples, &mixRSamples, 2);

			//Calculate distances between results and originals
			/*double dist = euclideanSquaredDistance(&(separated->data()[0][0]), &ilSamples);
			dist += euclideanSquaredDistance(&(separated->data()[0][1]), &irSamples);
			dist += euclideanSquaredDistance(&(separated->data()[1][0]), &jlSamples);
			dist += euclideanSquaredDistance(&(separated->data()[1][1]), &jrSamples);
			double distSwitched = euclideanSquaredDistance(&(separated->data()[0][0]), &jlSamples);
			distSwitched += euclideanSquaredDistance(&(separated->data()[0][1]), &jrSamples);
			distSwitched += euclideanSquaredDistance(&(separated->data()[1][0]), &ilSamples);
			distSwitched += euclideanSquaredDistance(&(separated->data()[1][1]), &irSamples);*/
			/*double dist = proportionalDifference(&(separated->data()[0][0]), &(separated->data()[0][1]), &ilSamples, &irSamples);
			dist += proportionalDifference(&(separated->data()[1][0]), &(separated->data()[1][1]), &jlSamples, &jrSamples);
			double distSwitched = proportionalDifference(&(separated->data()[0][0]), &(separated->data()[0][1]), &jlSamples, &jrSamples);
			distSwitched += proportionalDifference(&(separated->data()[1][0]), &(separated->data()[1][1]), &ilSamples, &irSamples);*/
			double dist = stereoCosineSimilarity(&(separated->data()[0][0]), &(separated->data()[0][1]), &ilSamples, &irSamples);
			dist += stereoCosineSimilarity(&(separated->data()[1][0]), &(separated->data()[1][1]), &jlSamples, &jrSamples);
			double distSwitched = stereoCosineSimilarity(&(separated->data()[0][0]), &(separated->data()[0][1]), &jlSamples, &jrSamples);
			distSwitched += stereoCosineSimilarity(&(separated->data()[1][0]), &(separated->data()[1][1]), &ilSamples, &irSamples);
			double bestDistance = ((dist>distSwitched)?dist:distSwitched)/2;
			sumDistance += bestDistance;

			cout << bestDistance << endl;

			delete separated;
		}
	}
	cout << endl;

	int numOfPairs = (sampleFiles.size() * (sampleFiles.size() - 1)) / 2;
	return sumDistance / numOfPairs;
}

void SinusoidalModelSeparation::optimiseParams(void)
{
	/*double mu = 4e-23;
	while(mu > 1e-24)
	{
		cout << "WA = " << SinusoidalTrajectory::distanceWeightAmp << endl;
		cout << "WH = " << SinusoidalTrajectory::distanceWeightHarm << endl;
		cout << "WS = " << SinusoidalTrajectory::distanceWeightStereo << endl;
		cout << "WP = " << SinusoidalTrajectory::distanceMissPenalty << endl;

		double currentError = totalErrorFromTestSet();

		cout << "Error = " << currentError << endl;

		SinusoidalTrajectory::distanceWeightAmp += 0.1;
		double newError = totalErrorFromTestSet();
		double derrdwa = (newError - currentError)/0.1;
		SinusoidalTrajectory::distanceWeightAmp -= 0.1;
		cout << "derr/dwa = " << derrdwa << endl;

		SinusoidalTrajectory::distanceWeightHarm += 0.1;
		newError = totalErrorFromTestSet();
		double derrdwh = (newError - currentError)/0.1;
		SinusoidalTrajectory::distanceWeightHarm -= 0.1;
		cout << "derr/dwh = " << derrdwh << endl;

		SinusoidalTrajectory::distanceWeightStereo += 0.1;
		newError = totalErrorFromTestSet();
		double derrdws = (newError - currentError)/0.1;
		SinusoidalTrajectory::distanceWeightStereo -= 0.1;
		cout << "derr/dws = " << derrdws << endl;

		SinusoidalTrajectory::distanceMissPenalty += 0.1;
		newError = totalErrorFromTestSet();
		double derrdwp = (newError - currentError)/0.1;
		SinusoidalTrajectory::distanceMissPenalty -= 0.1;
		cout << "derr/dwp = " << derrdwp << endl;
		*/
		/*SinusoidalTrajectory::distanceWeightAmp -= mu * derrdwa;
		SinusoidalTrajectory::distanceWeightHarm -= mu * derrdwh;
		SinusoidalTrajectory::distanceWeightStereo -= mu * derrdws;
		SinusoidalTrajectory::distanceMissPenalty -= mu * derrdwp;*/
	/*
		SinusoidalTrajectory::distanceWeightAmp -= .1 * derrdwa / abs(derrdwa);
		SinusoidalTrajectory::distanceWeightHarm -= .1 * derrdwh / abs(derrdwh);
		SinusoidalTrajectory::distanceWeightStereo -= .1 * derrdws / abs(derrdws);
		SinusoidalTrajectory::distanceMissPenalty -= .1 * derrdwp / abs(derrdwp);

		mu *= 0.9;
	}*/
	
	/*double minError = numeric_limits<double>::max();
	double bestWA = 0.;
	for(double wa=.4; wa<=1.2; wa += .1)
	{
		SinusoidalTrajectory::distanceWeightAmp = wa;
		cout << "WA = " << wa << "; ";
		double currentError = totalErrorFromStereoTestSet();
		cout << "Error = " << currentError << endl;
		if(minError > currentError)
		{
			minError = currentError;
			bestWA = wa;
		}
	}
	cout << "Best Error = " << minError << endl;
	cout << "Best WA = " << bestWA << endl;*/

	/*double minError = numeric_limits<double>::max();
	double bestWH = 0.;
	for(double wh=900; wh<=1000; wh += 10)
	{
		SinusoidalTrajectory::distanceWeightHarm = wh;
		cout << "WH = " << wh << "; ";
		double currentError = totalErrorFromStereoTestSet();
		cout << "Error = " << currentError << endl;
		if(minError > currentError)
		{
			minError = currentError;
			bestWH = wh;
		}
	}
	cout << "Best Error = " << minError << endl;
	cout << "Best WH = " << bestWH << endl;*/

	/*double minError = numeric_limits<double>::max();
	double bestMP = 0.;
	for(double mp=0; mp<=30; mp += 3)
	{
		SinusoidalTrajectory::distanceMissPenalty = mp;
		cout << "MP = " << mp << "; ";
		double currentError = totalErrorFromStereoTestSet();
		cout << "Error = " << currentError << endl;
		if(minError > currentError)
		{
			minError = currentError;
			bestMP = mp;
		}
	}
	cout << "Best Error = " << minError << endl;
	cout << "Best MP = " << bestMP << endl;*/

	/*double minError = numeric_limits<double>::max();
	double bestWF = 0.;
	for(double wf=0; wf<=100; wf += 10)
	{
		SinusoidalTrajectory::distanceWeightFreq = wf;
		cout << "WF = " << wf << "; ";
		double currentError = totalErrorFromStereoTestSet();
		cout << "Error = " << currentError << endl;
		if(minError > currentError)
		{
			minError = currentError;
			bestWF = wf;
		}
	}
	cout << "Best Error = " << minError << endl;
	cout << "Best WF = " << bestWF << endl;*/

	/*double minError = numeric_limits<double>::max();
	double bestWS = 0.;
	for(double ws=30; ws<=50; ws += 1)
	{
		SinusoidalTrajectory::distanceWeightStereo = ws;
		cout << "WS = " << ws << "; ";
		double currentError = totalErrorFromStereoTestSet();
		cout << "Error = " << currentError << endl;
		if(minError > currentError)
		{
			minError = currentError;
			bestWS = ws;
		}
	}
	cout << "Best Error = " << minError << endl;
	cout << "Best WS = " << bestWS << endl;*/

	/*double minError = numeric_limits<double>::max();
	double bestWO = 0.;
	for(double wo=1.; wo<=3; wo += .1)
	{
		SinusoidalTrajectory::distanceWeightOnset = wo;
		cout << "WO = " << wo << "; ";
		double currentError = totalErrorFromStereoTestSet();
		cout << "Error = " << currentError << endl;
		if(minError > currentError)
		{
			minError = currentError;
			bestWO = wo;
		}
	}
	cout << "Best Error = " << minError << endl;
	cout << "Best WO = " << bestWO << endl;*/

	/*double minError = numeric_limits<double>::max();
	double bestDT = 0.;
	//for(double dt=0; dt<=10; dt += 0.5)
	for(double dt=3.2; dt<=10; dt += 10)
	{
		clusterDistThreshold = dt;
		cout << "DT = " << dt << "; ";
		double currentError = totalErrorFromStereoTestSet();
		cout << "Error = " << currentError << endl;
		if(minError > currentError)
		{
			minError = currentError;
			bestDT = dt;
		}
	}
	cout << "Best Error = " << minError << endl;
	cout << "Best DT = " << bestDT << endl;*/

	double minError = numeric_limits<double>::max();
	double bestS = 0.;
	for(double s=2e-4; s<=3e-4; s += 1e-5)
	{
		softClusterStiffness = s;
		cout << "Stiffness = " << s << "; ";
		double currentError = totalErrorFromStereoTestSet();
		cout << "Error = " << currentError << endl;
		if(minError > currentError)
		{
			minError = currentError;
			bestS = s;
		}
	}
	cout << "Best Error = " << minError << endl;
	cout << "Best Stiffness = " << bestS << endl;
}
