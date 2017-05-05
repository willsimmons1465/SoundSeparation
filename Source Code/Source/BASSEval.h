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
	static vector<WavFileManager*>* testSetSounds();
	static BASSStats stereoTest(vector<vector<Eigen::VectorXd>>* originals, Separation::FeatureOption fOp, Separation::ClusterOption cOp, bool reversible);
	static void pairWiseStereoTest(vector<WavFileManager*>* sampleFiles, Separation::FeatureOption fOp, Separation::ClusterOption cOp, bool reversible, bool fullResults);
	static void standardTest(Separation::FeatureOption fOp, Separation::ClusterOption cOp, bool reversible);
	static void incrementKTest(Separation::FeatureOption fOp);
	static void incrementStereoTest(Separation::FeatureOption fOp);
	static void incrementOffsetTest(Separation::FeatureOption fOp);
	static void incrementFrequencyTest(Separation::FeatureOption fOp);
	static void incrementPitchTest(Separation::FeatureOption fOp);
	static void noiseTest(Separation::FeatureOption fOp, double noiseMult);
};
