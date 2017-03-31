/**
* CLASS: SinusoidalTrajectory
* AUTHOR: WVS22
*
* Instances of this class represent individual sinusoids with frequency, amplitude and stereo envelopes.
*/
#pragma once
#include "SinusoidalTrajectoryPoint.h"
#include <vector>
using namespace std;
class SinusoidalTrajectory
{
	static const double trajectoryDeltaFmax;
public:
	static double distanceWeightFreq; //was const value of 1 to give normalised weights
	static double distanceWeightAmp; //not const as can be modified in SinusoidalModelSeparation::optimiseParams
	static double distanceWeightHarm;
	static double distanceWeightStereo;
	static double distanceWeightOnset;
	static double distanceMissPenalty;
	static double distanceMissPenaltyFreq;
	static double distanceMissPenaltyAmp;
	static double distanceMissPenaltyStereo;
	vector<SinusoidalTrajectoryPoint> data;
	long startIndex;
	long endIndex;
	double meanFreq;
	double meanAmp;
	SinusoidalTrajectory(long start, SinusoidalTrajectoryPoint firstPoint) : startIndex(start), endIndex(start), data() {data.push_back(firstPoint);}
	
	/**
	* SinusoidalTrajectory::addPoint(SinusoidalTrajectoryPoint)
	*
	* Adds a new point to the end of the trajectory.
	* Assumes each call to this method should be made for the time frame immediately
	* following the previous call, so all points are consecutive in time and once it
	* has ended, no more points are added.
	*/
	inline void addPoint(SinusoidalTrajectoryPoint point){data.push_back(point); endIndex++;}
	
	/**
	* SinusoidalTrajectory::canAccept(SinusoidalTrajectoryPoint)
	* 
	* Returns true iff point could be added to the end of the trajectory without
	* a significant jump in frequency. trajectoryDeltaFmax is the maximum percentage
	* change in frequency between adjacent time frames.
	*/
	bool canAccept(SinusoidalTrajectoryPoint point);

	/**
	* SinusoidalTrajectory::normalise()
	*
	* Normalises the amp and freq values of the points on the trajectory to have
	* means of 1. Does not affect the stereo position or trueFreq fields, so original
	* frequency can be read from trueFreq.
	* This should be called exactly once after all points have been added to the trajectory,
	* otherwise the meanFreq field may not be correct.
	*/
	void normalise(void);

	/**
	* SinusoidalTrajectory::distance(SinusoidalTrajectory, SinusoidalTrajectory)
	*
	* Determines the "distance" between two trajectories based on Euclidean distance
	* between the amplitude, frequency and stereo vectors over the overlapping time period.
	* A measure of harmonic distance is added according to Virtanen's equation.
	* A penalty figure is added instead of amplitude, frequency and stereo distance if
	* trajectories do not overlap in time.
	*/
	static double distance(SinusoidalTrajectory sin0, SinusoidalTrajectory sin1, double minFreq);

	static double distanceMiss(SinusoidalTrajectory sin0, SinusoidalTrajectory sin1);
	static double distanceHarm(SinusoidalTrajectory sin0, SinusoidalTrajectory sin1, double minFreq);
	static double distanceStereo(SinusoidalTrajectory sin0, SinusoidalTrajectory sin1);
	static double distanceOnset(SinusoidalTrajectory sin0, SinusoidalTrajectory sin1);
	static double distanceAmp(SinusoidalTrajectory sin0, SinusoidalTrajectory sin1);
	static double distanceFreq(SinusoidalTrajectory sin0, SinusoidalTrajectory sin1);
};
