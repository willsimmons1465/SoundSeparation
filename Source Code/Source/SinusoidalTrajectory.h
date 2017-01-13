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
	static const double distanceWeightFreq;
	static const double distanceWeightAmp;
	static const double distanceWeightHarm;
	static const double distanceWeightStereo;
	static const double distanceMissPenalty;
public:
	vector<SinusoidalTrajectoryPoint> data;
	long startIndex;
	long endIndex;
	double meanFreq;
	SinusoidalTrajectory(long start, SinusoidalTrajectoryPoint firstPoint) : startIndex(start), endIndex(start), data() {data.push_back(firstPoint);}
	inline void addPoint(SinusoidalTrajectoryPoint point){data.push_back(point); endIndex++;}
	bool canAccept(SinusoidalTrajectoryPoint point);
	void normalise();
	static double distance(SinusoidalTrajectory sin0, SinusoidalTrajectory sin1, double minFreq);

};
