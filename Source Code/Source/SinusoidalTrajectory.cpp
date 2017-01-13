#include "SinusoidalTrajectory.h"

const double SinusoidalTrajectory::trajectoryDeltaFmax = 0.02;
const double SinusoidalTrajectory::distanceWeightFreq = 1.;
const double SinusoidalTrajectory::distanceWeightAmp = 2.5;
const double SinusoidalTrajectory::distanceWeightHarm = 1000.;
const double SinusoidalTrajectory::distanceWeightStereo = 1.;
const double SinusoidalTrajectory::distanceMissPenalty = 1.;

bool SinusoidalTrajectory::canAccept(SinusoidalTrajectoryPoint point)
{
	SinusoidalTrajectoryPoint lastPoint = data.back();
	return abs(log(point.freq/lastPoint.freq)) < log(1+trajectoryDeltaFmax);
}

void SinusoidalTrajectory::normalise(void)
{
	//Calculate mean amp and freq
	double sumAmp = 0.;
	double sumFreq = 0.;
	for(long i=0; i<data.size(); i++)
	{
		sumAmp += data[i].amp;
		sumFreq += data[i].freq;
	}
	const double meanAmp = sumAmp/data.size();
	meanFreq = sumFreq/data.size();

	//Scale all points
	for(long i=0; i<data.size(); i++)
	{
		data[i].amp /= meanAmp;
		data[i].freq /= meanFreq;
	}
}


double SinusoidalTrajectory::distance(SinusoidalTrajectory sin0, SinusoidalTrajectory sin1, double minFreq)
{
	double df = 0.;
	double da = 0.;
	double ds = 0.;
	bool overlap = true;

	//Identify overlapping region
	long t0;
	if(sin0.startIndex < sin1.startIndex)
	{
		t0 = sin1.startIndex;
	}
	else
	{
		t0 = sin0.startIndex;
	}
	long t1;
	if(sin0.endIndex < sin1.endIndex)
	{
		t1 = sin0.endIndex;
	}
	else
	{
		t1 = sin1.endIndex;
	}

	//Handle overlap penalty/calculate amp, freq + stereo distances
	if(t1 < t0)
	{
		overlap = false;
	}
	else
	{
		for(long t=t0; t<=t1; t++)
		{
			df += pow(sin0.data[t - sin0.startIndex].freq - sin1.data[t - sin1.startIndex].freq, 2);
			da += pow(sin0.data[t - sin0.startIndex].amp - sin1.data[t - sin1.startIndex].amp, 2);
			ds += pow(sin0.data[t - sin0.startIndex].stereo - sin1.data[t - sin1.startIndex].stereo, 2);
		}
		double scale = 1/(1+t1-t0);
		df *= scale;
		da *= scale;
		ds *= scale;
	}

	//Calculate the harmonic distance
	//Used ceil instead of floor to cover error case
	//e.g. let minFreq = sin0.meanFreq = 5 and sin1.meanFreq = 24.9
	//We want to consider a = 1, b = 5
	long maxA = ceil(sin0.meanFreq/minFreq);
	long maxB = ceil(sin1.meanFreq/minFreq);
	double dh = numeric_limits<double>::max();
	for(long a = 1; a <= maxA; a++)
	{
		//For each value of a, the minimum values will be from one of the two b
		//values around expectedB
		double expectedB = a*sin1.meanFreq/sin0.meanFreq;

		//expectedB increases with a, so at this point there can be no more options
		if(expectedB > maxB)
		{
			break;
		}

		//Calculate the distance for each value of b
		double dhab = abs(log(sin0.meanFreq*floor(expectedB)/(sin1.meanFreq*a)));
		if(dh > dhab)
		{
			dh = dhab;
		}
		dhab = abs(log(sin0.meanFreq*ceil(expectedB)/(sin1.meanFreq*a)));
		if(dh > dhab)
		{
			dh = dhab;
		}
	}

	//Return summed distance
	return distanceWeightFreq*df + distanceWeightAmp*da + distanceWeightHarm*dh + distanceWeightStereo*ds + (overlap?0.:distanceMissPenalty);
}
