#include "Separation.h"
#include "Cluster.h"
#include "SinusoidalTrajectory.h"
#include "SinusoidalTrajectoryPoint.h"
#include "Transform.h"
#include <complex>
#include <iostream>

//const double Separation::noiseThreshold = 1e20;
//const double Separation::noiseThreshold = 2.5e21;
//const double Separation::noiseThreshold = 1e136;
//const double Separation::noiseThreshold = 1e110;
//double Separation::softClusterStiffness = 0.00022;
const double naiveClusterThreshold = 3.2;
//double Separation::matrixWeightSource = 1.;
//double Separation::matrixWeightMix = 1.;

//const double Separation::noiseThreshold = 2.5e21;
const double Separation::noiseThreshold = 0.05e21;
double Separation::softClusterStiffnessSin = 0.0001;
double Separation::softClusterStiffnessMat = 15;
double Separation::matrixWeightMix = 0.002;
double Separation::matrixWeightSource = 0.5;

inline double quadraticFitTurningPoint(double y_1, double y0, double y1)
{
	return -0.5*(y1-y_1)/(y1+y_1-2*y0);
}

vector<double>* peakPositions(Eigen::VectorXd* values)
{
	vector<double>* retVal = new vector<double>();
	double prev = 0.;
	double current = (*values)(0);
	double next = (*values)(1);
	for(long i=2; i<values->size(); i++)
	{
		prev = current;
		current = next;
		next = (*values)(i, 0);
		if(current > prev && current > next)
		{
			retVal->push_back(i-1+quadraticFitTurningPoint(prev, current, next));
		}
	}

	return retVal;
}

double approxPowerAroundPosition(Eigen::VectorXd* values, double position)
{
	long lPos = (long) position;
	//double part = position-lPos;
	//return (1-part)*(*values)(lPos)*(*values)(lPos) + part*(*values)(lPos+1)*(*values)(lPos+1);
	return (*values)(lPos) + (*values)(lPos+1);
}

double approxPowerAroundPosition(Eigen::VectorXcd* values, double position)
{
	long lPos = (long) position;
	//double part = position-lPos;
	//return (1-part)*norm((*values)(lPos)) + part*norm((*values)(lPos+1));
	return norm((*values)(lPos)) + norm((*values)(lPos+1));
}

double stereoScaleAroundPosition(Eigen::VectorXcd* left, Eigen::VectorXcd* right, double position)
{
	double lPower = approxPowerAroundPosition(left, position);
	double rPower = approxPowerAroundPosition(right, position);
	return lPower / (lPower + rPower);
}

vector<SinusoidalTrajectory>* extractSinusoids(Eigen::MatrixXd* absSpec, Eigen::MatrixXcd* lSpectrum, Eigen::MatrixXcd* rSpectrum)
{
	long tCount = absSpec->rows();
	vector<SinusoidalTrajectory>* sins = new vector<SinusoidalTrajectory>();
	for(long t=0; t<tCount; t++)
	{
		//Find local peak frequencies for current time frame
		Eigen::VectorXd frame = absSpec->row(t);
		Eigen::VectorXcd lframe = lSpectrum->row(t);
		Eigen::VectorXcd rframe = rSpectrum->row(t);
		vector<double>* peakPos = peakPositions(&frame);
		for(long i=0; i<peakPos->size(); i++)
		{
			//Create a trajectory point for each peak
			double p = (*peakPos)[i];
			SinusoidalTrajectoryPoint point(approxPowerAroundPosition(&frame,p),
				p,
				stereoScaleAroundPosition(&lframe, &rframe, p));
			
			//Extend an existing trajectory if possible
			bool newSin = true;
			for(long i=0; i<sins->size(); i++)
			{
				if((*sins)[i].endIndex == t-1 && (*sins)[i].canAccept(point))
				{
					(*sins)[i].addPoint(point);
					newSin = false;
					break;
				}
			}

			//If no trajectory could be extended, create a new one
			if(newSin)
			{
				sins->push_back(SinusoidalTrajectory(t,point));
			}
		}

		delete peakPos;
	}

	return sins;
}

Eigen::MatrixXd* distanceMatrixFromSinusoids(vector<SinusoidalTrajectory>* sins)
{
	long numOfSins = sins->size();
	double minFreq = numeric_limits<double>::max();
	for(long i=0; i<numOfSins; i++)
	{
		(*sins)[i].normalise();
		if(minFreq > (*sins)[i].meanFreq)
		{
			minFreq = (*sins)[i].meanFreq;
		}
	}

	Eigen::MatrixXd* d = new Eigen::MatrixXd(numOfSins, numOfSins*6);
	for(long i=0; i<numOfSins; i++)
	{
		for(long j=i; j<numOfSins; j++)
		{
			(*d)(i, 6*j) = SinusoidalTrajectory::distanceHarm((*sins)[i], (*sins)[j], minFreq);
			(*d)(i, 6*j+1) = SinusoidalTrajectory::distanceStereo((*sins)[i], (*sins)[j]);
			(*d)(i, 6*j+2) = SinusoidalTrajectory::distanceOnset((*sins)[i], (*sins)[j]);
			(*d)(i, 6*j+3) = SinusoidalTrajectory::distanceMiss((*sins)[i], (*sins)[j]);
			(*d)(i, 6*j+4) = SinusoidalTrajectory::distanceAmp((*sins)[i], (*sins)[j]);
			(*d)(i, 6*j+5) = SinusoidalTrajectory::distanceFreq((*sins)[i], (*sins)[j]);

			if(j!=i)
			{
				(*d)(j, 6*i) = (*d)(i, 6*j);
				(*d)(j, 6*i+1) = (*d)(i, 6*j+1);
				(*d)(j, 6*i+2) = (*d)(i, 6*j+2);
				(*d)(j, 6*i+3) = (*d)(i, 6*j+3);
				(*d)(j, 6*i+4) = (*d)(i, 6*j+4);
				(*d)(j, 6*i+5) = (*d)(i, 6*j+5);
			}
		}
	}

	return d;
}

Eigen::MatrixXd* distanceMatrixFromMatrixFactors(Eigen::MatrixXd* mix, Eigen::MatrixXd* source)
{
	long factorCount = mix->cols();
	Eigen::MatrixXd* features = new Eigen::MatrixXd(factorCount, factorCount);
	for(int i=0; i<factorCount; i++)
	{
		for(int j=0; j<i; j++)
		{
			(*features)(i, j) = (*features)(j, i);
			(*features)(i, factorCount+j) = (*features)(j, factorCount+i);
		}
		for(int j=i; j<factorCount; j++)
		{
			(*features)(i, j) = (mix->col(i) - mix->col(j)).norm();
			(*features)(i, factorCount+j) = (source->row(i) - source->row(j)).norm();
		}
	}

	return features;
}

Eigen::MatrixXd* naiveSinusoidClustering(vector<SinusoidalTrajectory>* sins, int numOfSources)
{
	long featureCount = sins->size();
	double minFreq = numeric_limits<double>::max();
	for(long i=0; i<featureCount; i++)
	{
		if(minFreq > (*sins)[i].meanFreq)
		{
			minFreq = (*sins)[i].meanFreq;
		}
	}
	vector<int> clusterTags = vector<int>(featureCount);
	for(long i=0; i<featureCount; i++)
	{
		clusterTags[i] = -1;
	}
	for(int s=0; s<numOfSources; s++)
	{
		long seed = -1;
		double maxAmp = 0.;
		for(long i=0; i<featureCount; i++)
		{
			if(clusterTags[i] == -1 && maxAmp < (*sins)[i].meanAmp)
			{
				seed = i;
				maxAmp = (*sins)[i].meanAmp;
			}
		}
		if (seed == -1)
		{
			break;
		}
		for(long i=0; i<featureCount; i++)
		{
			if(SinusoidalTrajectory::distanceHarm((*sins)[i], (*sins)[seed], minFreq) < naiveClusterThreshold)
			{
				clusterTags[i] = s;
			}
		}
	}
	for(long i=0; i<featureCount; i++)
	{
		if(clusterTags[i] != -1)
		{
			continue;
		}
		vector<double> maxHarmDists = vector<double>(numOfSources);
		for(long j=0; j<featureCount; j++)
		{
			int jTag = clusterTags[j];
			if(jTag != -1)
			{
				double dh = SinusoidalTrajectory::distanceHarm((*sins)[i], (*sins)[j], minFreq);
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
				clusterTags[i] = s;
			}
		}
	}

	Eigen::MatrixXd* retVal = new Eigen::MatrixXd(featureCount, numOfSources);
	(*retVal) = Eigen::MatrixXd::Zero(featureCount, numOfSources);
	for(long i=0; i<featureCount; i++)
	{
		(*retVal)(i, clusterTags[i]) = 1.;
	}

	return retVal;
}

vector<vector<Eigen::VectorXd>>* Separation::separate(Eigen::VectorXd* lSampleVector, Eigen::VectorXd* rSampleVector,
	int numOfSources, FeatureOption fOp, ClusterOption cOp, bool reversible, bool verbose)
{
	//Obtain spectra
	Eigen::MatrixXcd* lSpectrum = Transform::stft(lSampleVector, windowSize, hopSize, nfft);
	Eigen::MatrixXcd* rSpectrum = Transform::stft(rSampleVector, windowSize, hopSize, nfft);

	if (verbose)
	{
		cout << "Spectra obtained" << endl;
	}

	//Get sizes
	long tCount = lSpectrum->rows();
	long fCount = lSpectrum->cols();

	//Get overall power spectrum
	Eigen::MatrixXd absSpec = (lSpectrum->cwiseAbs2() + rSpectrum->cwiseAbs2());

	if (verbose)
	{
		cout << "Power spectrum obtained" << endl;
	}

	//Threshold power spectrum
	if (fOp == SINUSOIDS)
	{
		for(long t=0; t<tCount; t++)
		{
			for(long f=0; f<fCount; f++)
			{
				if(noiseThreshold > absSpec(t, f))
				{
					absSpec(t, f) = 0;
				}
			}
		}

		if (verbose)
		{
			cout << "Power spectrum thresholded" << endl;
		}
	}

	//Obtain features
	vector<SinusoidalTrajectory>* sins;
	Eigen::MatrixXd mix;
	Eigen::MatrixXd source;
	switch (fOp)
	{
	case SINUSOIDS:
		sins = extractSinusoids(&absSpec, lSpectrum, rSpectrum);

		if (verbose)
		{
			cout << sins->size() << " Sinusoids extracted" << endl;
		}

		cout << sins->size() << endl;

		break;
	case MATRIXFACTORS:
		vector<Eigen::MatrixXd>* nmfResults = Transform::nmf(&absSpec, nmfFactors);
		mix = nmfResults->data()[0];
		source = nmfResults->data()[1];
		delete nmfResults;

		if (verbose)
		{
			cout << "NMF computed" << endl;
		}
	}

	//Construct feature vectors for clustering
	Eigen::MatrixXd* featureVectors;
	switch (fOp)
	{
	case SINUSOIDS:
		featureVectors = distanceMatrixFromSinusoids(sins);

		if (verbose)
		{
			cout << "Distance matrix calculated" << endl;
		}
		break;
	case MATRIXFACTORS:
		//featureVectors = distanceMatrixFromMatrixFactors(&mix, &source);
		featureVectors = new Eigen::MatrixXd(nmfFactors, tCount + fCount);
		for(int i=0; i<nmfFactors; i++)
		{
			featureVectors->row(i) << (source.row(i).normalized() * matrixWeightSource), (mix.transpose().row(i).normalized() * matrixWeightMix);
		}
		//(*featureVectors) << (source * matrixWeightSource), (mix.transpose() * matrixWeightMix);

		if (verbose)
		{
			cout << "Feature vectors created" << endl;
		}
	}

	//Cluster features
	Eigen::MatrixXd* clusterWeights;
	long featureCount = featureVectors->rows();
	switch (cOp)
	{
	case HARD:
		{
		vector<int>* clusterTags = Cluster::kClusterLloyd(featureVectors, numOfSources);
		clusterWeights = new Eigen::MatrixXd(featureCount, numOfSources);
		(*clusterWeights) = Eigen::MatrixXd::Zero(featureCount, numOfSources);
		for(long i=0; i<featureCount; i++)
		{
			(*clusterWeights)(i, (*clusterTags)[i]) = 1.;
		}
		delete clusterTags;
		}

		break;
	case SOFT:
		switch (fOp)
		{
		case SINUSOIDS:
			clusterWeights = Cluster::softKCluster(featureVectors, numOfSources, softClusterStiffnessSin);
			break;
		case MATRIXFACTORS:
			clusterWeights = Cluster::softKCluster(featureVectors, numOfSources, softClusterStiffnessMat);
		}

		break;
	case MATRIX:
		{
		vector<Eigen::MatrixXd>* nmfResults = Transform::nmf(featureVectors, numOfSources);
		clusterWeights = new Eigen::MatrixXd(featureCount, numOfSources);
		(*clusterWeights) = (*nmfResults)[0];
		for(long i=0; i<featureCount; i++)
		{
			if(clusterWeights->row(i).isZero())
			{
				clusterWeights->row(i).setOnes();
			}
			clusterWeights->row(i).normalize();
		}
		delete nmfResults;
		}
		
		break;
	case NAIVE:
		clusterWeights = naiveSinusoidClustering(sins, numOfSources);
	}

	delete featureVectors;

	if (verbose)
	{
		cout << "Features clustered" << endl;
	}

	//Construct each sound's spectrum
	vector<Eigen::MatrixXcd> lSoundSpecs = vector<Eigen::MatrixXcd>();
	vector<Eigen::MatrixXcd> rSoundSpecs = vector<Eigen::MatrixXcd>();
	switch (fOp)
	{
	case SINUSOIDS:
		//cout << (*clusterWeights) << endl;
		for(int s=0; s<numOfSources; s++)
		{
			//Initialise spectra as zero-matrices
			lSoundSpecs.push_back(Eigen::MatrixXcd::Zero(tCount, fCount));
			rSoundSpecs.push_back(Eigen::MatrixXcd::Zero(tCount, fCount));

			//Add peaks from trajectories in the current cluster
			for(long i=0; i<featureCount; i++)
			{
				for(long t=(*sins)[i].startIndex; t<=(*sins)[i].endIndex; t++)
				{
					long freq = (long) (*sins)[i].data[t-(*sins)[i].startIndex].trueFreq;
					lSoundSpecs[s](t, freq) = (*lSpectrum)(t, freq) * (*clusterWeights)(i, s);
					lSoundSpecs[s](t, freq+1) = (*lSpectrum)(t, freq+1) * (*clusterWeights)(i, s);
					rSoundSpecs[s](t, freq) = (*rSpectrum)(t, freq) * (*clusterWeights)(i, s);
					rSoundSpecs[s](t, freq+1) = (*rSpectrum)(t, freq+1) * (*clusterWeights)(i, s);
				}
			}
		}
		break;
	case MATRIXFACTORS:
		for(int s=0; s<numOfSources; s++)
		{
			//Obtain absolute spectrum for source
			Eigen::MatrixXd soundMix = mix;
			Eigen::MatrixXd soundSource = source;
			for(int i=0; i<featureCount; i++)
			{
				soundMix.col(i) *= (*clusterWeights)(i, s);
				soundSource.row(i) *= (*clusterWeights)(i, s);
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

			lSoundSpecs.push_back(lSoundSpec);
			rSoundSpecs.push_back(rSoundSpec);
		}

	}

	if (fOp == SINUSOIDS)
	{
		delete sins;
	}

	if(verbose)
	{
		cout << "Reconstructed spectra" << endl;
	}

	if(reversible)
	{
		for(int s=0; s<numOfSources; s++)
		{
			(*lSpectrum) -= lSoundSpecs[s];
			(*rSpectrum) -= rSoundSpecs[s];
		}
		for(int s=0; s<numOfSources; s++)
		{
			lSoundSpecs[s] += (*lSpectrum) / (double) numOfSources;
			rSoundSpecs[s] += (*rSpectrum) / (double) numOfSources;
		}

		if(verbose)
		{
			cout << "Added non-feature sounds" << endl;
		}
	}

	//Reconstruct audio
	vector<vector<Eigen::VectorXd>>* retVal = new vector<vector<Eigen::VectorXd>>();
	for(int s=0; s<numOfSources; s++)
	{
		//Calculate the samples from the spectra
		Eigen::VectorXd* lSamples = Transform::istft(&(lSoundSpecs[s]), hopSize);
		Eigen::VectorXd* rSamples = Transform::istft(&(rSoundSpecs[s]), hopSize);
		retVal->push_back(vector<Eigen::VectorXd>());
		retVal->back().push_back(*lSamples);
		retVal->back().push_back(*rSamples);

		delete lSamples;
		delete rSamples;
	}

	if (verbose)
	{
		cout << "Reconstructed audio" << endl;
	}

	//Clean up
	delete lSpectrum;
	delete rSpectrum;
	delete clusterWeights;

	return retVal;
}
