/**
* CLASS: BASSEval
* AUTHOR: WVS22
*
* Defines a collection of static methods for evaluating the performance of the separation algorithm
*/
#pragma once
#include <vector>
#include <Eigen/Core>
#include "WavFileManager.h"
#include "Separation.h"
using namespace std;
class BASSEval
{
public:
	struct BASSStats
	{
		double avgSDR;
		double avgSIR;
		double avgSAR;
	};
	static vector<WavFileManager*>* validationSetSounds();
	static BASSStats stereoTest(vector<vector<Eigen::VectorXd>>* originals, Separation::FeatureOption fOp, Separation::ClusterOption cOp, bool reversible);
	static void pairWiseStereoTest(vector<WavFileManager*>* sampleFiles, Separation::FeatureOption fOp, Separation::ClusterOption cOp, bool reversible, bool fullResults);
};
