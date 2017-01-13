#include "Cluster.h"
//#include <cstdlib>

//Helper method for kClusterLloyd
//Returns square of Euclidean distance between vectors.
double distance(vector<double>* vector0, vector<double>* vector1)
{
	double sum = 0.;
	for(long f=0; f<vector0->size(); f++)
	{
		sum += pow((*vector0)[f] - (*vector1)[f], 2);
	}
	return sum;
}

vector<int>* Cluster::kClusterLloyd(vector<vector<double>>* featureVectors, int k)
{
	//Get parameters
	long vectorCount = featureVectors->size();
	long featureCount = (*featureVectors)[0].size();

	//Identify ranges of values for each feature
	/*vector<double> mins = vector<double>();
	vector<double> maxs = vector<double>();
	for(long f=0; f<featureCount; f++)
	{
		mins.push_back((*featureVectors)[0][f]);
		maxs.push_back((*featureVectors)[0][f]);
	}
	for(long v=1; v<vectorCount; v++)
	{
		for(long f=0; f<featureCount; f++)
		{
			double val = (*featureVectors)[v][f];
			if(val < mins[f])
			{
				mins[f] = val;
			}
			else if(val > maxs[f])
			{
				maxs[f] = val;
			}
		}
	}*/

	//Select random points as initial centres
	vector<vector<double>> centres = vector<vector<double>>();
	for(int c=0; c<k; c++)
	{
		long centreChoice = rand() % vectorCount;
		centres.push_back(vector<double>((*featureVectors)[centreChoice]));
	}

	//Improve centres
	vector<int>* clusterTags = new vector<int>(vectorCount);
	long tagsChanged = lloydConvergeLimit + 1;
	while (tagsChanged > lloydConvergeLimit)
	{
		//Cluster points to nearest centre
		tagsChanged = 0;
		for(long v=0; v<vectorCount; v++)
		{
			int bestTag = -1;
			double bestDistance = numeric_limits<double>::max();
			for(int c=0; c<k; c++)
			{
				double dist = distance(&((*featureVectors)[v]), &(centres[c]));
				if(bestDistance > dist)
				{
					bestDistance = dist;
					bestTag = c;
				}
			}
			if(bestTag != (*clusterTags)[v])
			{
				(*clusterTags)[v] = bestTag;
				tagsChanged++;
			}
		}
		
		//Recalculate the centres given their clusters
		for(int c=0; c<k; c++)
		{
			long clusterSize = 0;
			for(long f=0; f<featureCount; f++)
			{
				centres[c][f] = 0.;
			}
			for(long v=0; v<vectorCount; v++)
			{
				if(c == (*clusterTags)[v])
				{
					clusterSize++;
					for(long f=0; f<featureCount; f++)
					{
						centres[c][f] += (*featureVectors)[v][f];
					}
				}
			}
			for(long f=0; f<featureCount; f++)
			{
				centres[c][f] /= clusterSize;
			}
		}
	}

	return clusterTags;
}
