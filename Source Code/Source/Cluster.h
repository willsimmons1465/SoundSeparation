/**
* CLASS: Cluster
* AUTHOR: WVS22
*
* Defines a collection of static methods for clustering data.
*/
#pragma once
#include <vector>
using namespace std;
class Cluster
{
	static const long lloydConvergeLimit = 1;
public:
	static vector<int>* kClusterLloyd(vector<vector<double>>* featureVectors, int k);
};
