/**
* CLASS: SinusoidalTrajectoryPoint
* AUTHOR: WVS22
*
* Instances of this class represent individual points for a SinusoidalTrejectory instance.
*/
#pragma once
struct SinusoidalTrajectoryPoint
{
	double amp;
	double freq;
	double trueFreq;
	double stereo;
	SinusoidalTrajectoryPoint(double iamp, double ifreq, double istereo) : amp(iamp), freq(ifreq), trueFreq(ifreq), stereo(istereo) {}
};
