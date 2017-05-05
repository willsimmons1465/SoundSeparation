/**
* CLASS: Separation
* AUTHOR: WVS22
*
* Creates a static access point to the sound separation options.
*/
#pragma once
#include <vector>
#include <Eigen/Core>
using namespace std;
class Separation
{
	static const long windowSize = 8192;
	static const long hopSize = 1024;
	static const long nfft = 8192;
	static const long nmfFactors = 10;
	//static const long nmfFactors = 221;
	static const double noiseThreshold;
public:
	static double softClusterStiffnessSin;
	static double softClusterStiffnessMat;
	static double matrixWeightSource;
	static double matrixWeightMix;
	enum FeatureOption {
		SINUSOIDS,
		MATRIXFACTORS
	};
	enum ClusterOption {
		HARD,
		SOFT,
		MATRIX,
		NAIVE
	};

	static vector<vector<Eigen::VectorXd>>* separate(Eigen::VectorXd* lSampleVector, Eigen::VectorXd* rSampleVector,
		int numOfSources, FeatureOption fOp, ClusterOption cOp, bool reversible, bool verbose);
};
