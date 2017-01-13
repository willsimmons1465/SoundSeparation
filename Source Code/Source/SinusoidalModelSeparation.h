/**
* CLASS: SinusoidalModelSeparation
* AUTHOR: WVS22
*
* Creates a static access point to the Sinusoidal Modelling method of sound separation.
*/
#pragma once
#include <vector>
using namespace std;
class SinusoidalModelSeparation
{
	static const long windowSize = 4096;
	static const long hopSize = 1024;
	static const long nfft = 4096;
public:

	/**
	* SinusoidalModelSepatation::separate(vector<double>*, vector<double>*)
	*
	* Analyses the given stereo sound to identify the most significant sinusoids, then
	* clusters these to identify reasonable approximations to the separate sound sources.
	* Callers must provide the number of sound sources present.
	* The return value is a vector of numOfSources many sources. Each source is a pair of
	* vectors of samples.
	*/
	static vector<vector<vector<double>>>* separate(vector<double>* lSampleVector, vector<double>* rSampleVector, int numOfSources);

};
