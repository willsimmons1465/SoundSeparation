/**
* CLASS: Cluster
* AUTHOR: WVS22
*
* Defines a collection of static methods for clustering data.
*/
#pragma once
#include <vector>
#include "Matrix.h"
using namespace std;
class Cluster
{
	static const long lloydConvergeLimit = 1;
public:

	/**
	* Cluster::kClusterLloyd(Matrix<double>*, int)
	*
	* Implementation of Lloyd's algorithm for clustering.
	* featureVectors is expected to be a vector of points, each point described by a
	* vector of features. k is the number of clusters to group the points into.
	* Determines distances by Euclidean distance and calculates centres using centre of mass.
	* NOTE: no checks are made to ensure that the initial set of centres does not contain
	* duplicates. This may cause some clusters to be empty. If this happens, use srand to
	* set a different seed for the random number generator.
	* Convergence is determined when the number of points which switch clusters drops
	* below lloydConvergeLimit. Raise this value to improve speed, sacrificing quality.
	*/
	static vector<int>* kClusterLloyd(vector<vector<double>>* featureVectors, int k);

	static Matrix<double>* softKCluster(vector<vector<double>>* featureVectors, int k, double stiffness);
};
