/**
* CLASS: NMFSeparation
* AUTHOR: WVS22
*
* Creates a static access point to the NMF method of sound separation.
*/
#pragma once
#include <vector>
using namespace std;
class NMFSeparation
{
	static const long windowSize = 4096;
	static const long hopSize = 1024;
	static const long nfft = 4096;
public:

	/**
	* NMFSepatation::separate(vector<double>*, vector<double>*)
	*
	* Performs Non-negative Matrix Factorisation on the spectrogram of the provided audio.
	* Inverse-transforms the component matrices, interpreting them as separate sources.
	* The return value is a vector of numOfSources many sources. Each source is a pair of
	* vectors of samples.
	*/
	static vector<vector<vector<double>>>* separate(vector<double>* lSampleVector, vector<double>* rSampleVector, int numOfSources);
};