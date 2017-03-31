#include "Cluster.h"
#include <math.h>
//#include <cstdlib>
//#include <iostream>

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

vector<int>* Cluster::kClusterLloyd(Eigen::MatrixXd* featureVectors, int k)
{
	//Get parameters
	long vectorCount = featureVectors->rows();
	long featureCount = featureVectors->cols();

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
	/*vector<vector<double>> centres = vector<vector<double>>();
	for(int c=0; c<k; c++)
	{
		long centreChoice = rand() % vectorCount;
		centres.push_back(vector<double>((*featureVectors)[centreChoice]));
	}*/
	Eigen::MatrixXd centres = Eigen::MatrixXd(k, featureCount);
	for(int c=0; c<k; c++)
	{
		long centreChoice = rand() % vectorCount;
		centres.row(c) = featureVectors->row(centreChoice);
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
			//Eigen::VectorXd currentVector = featureVectors->row(v);
			int bestTag = -1;
			double bestDistance = numeric_limits<double>::max();
			for(int c=0; c<k; c++)
			{
				//double dist = distance(&((*featureVectors)[v]), &(centres[c]));
				//Eigen::VectorXd currentCentre = centres.row(c);
				//double dist = (currentVector - currentCentre).squaredNorm();
				double dist = (featureVectors->row(v) - centres.row(c)).squaredNorm();
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
		/*for(int c=0; c<k; c++)
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
		}*/

		centres.setZero(centres.rows(), centres.cols());
		vector<int> clusterSizes = vector<int>(k);
		for(long v=0; v<vectorCount; v++)
		{
			int cluster = (*clusterTags)[v];
			centres.row(cluster) += (*featureVectors).row(v);
			clusterSizes[cluster]++;
		}
		for(int c=0; c<k; c++)
		{
			centres.row(c) /= clusterSizes[c];
		}
	}

	return clusterTags;
}


Eigen::MatrixXd* Cluster::softKCluster(Eigen::MatrixXd* featureVectors, int k, double stiffness)
{
	long vectorCount = featureVectors->rows();
	long featureCount = featureVectors->cols();
	Eigen::MatrixXd* forceMatrix = new Eigen::MatrixXd(vectorCount, k);

	Eigen::MatrixXd centres = Eigen::MatrixXd(k, featureCount);
	/*for(int c=0; c<k; c++)
	{
		long centreChoice = rand() % vectorCount;
		centres.push_back(vector<double>((*featureVectors)[centreChoice]));
	}*/
	for(int c=0; c<k; c++)
	{
		long centreChoice = rand() % vectorCount;
		centres.row(c) = featureVectors->row(centreChoice);
	}

	double lastError = numeric_limits<double>::max();
	while(true)
	{
		/*for(int c=0; c<k; c++)
		{
			cout << "Centre " << c << ": ";
			for(int f=0; f<featureCount; f++)
			{
				cout << centres[c][f] << ", ";
			}
			cout << endl;
		}*/

		double error = 0;
		//Eigen::VectorXd totalForCentre = Eigen::VectorXd(k);
		for(long v=0; v<vectorCount; v++)
		{
			//Eigen::VectorXd currentVector = featureVectors->row(v);
			//double totalForPoint = 0.;
			for(int c=0; c<k; c++)
			{
				//double dist = distance(&((*featureVectors)[v]), &(centres[c]));
				double dist = (featureVectors->row(v) - centres.row(c)).squaredNorm();
				(*forceMatrix)(v, c) = exp(-stiffness*dist);
				error += (*forceMatrix)(v, c) * dist;
				//totalForPoint += (*forceMatrix)(v, c);
			}
			/*for(int c=0; c<k; c++)
			{
				(*forceMatrix)(v, c) /= totalForPoint;
				totalForCentre[c] += (*forceMatrix)(v, c);
			}*/
			forceMatrix->row(v).normalize();
			
		}

		/*for(long v=0; v<vectorCount; v++)
		{
			cout << "Force on " << v << ": ";
			for(int c=0; c<k; c++)
			{
				cout << forceMatrix->data[v][c] << ", ";
			}
			cout << endl;
		}*/

		//if(abs(lastError-error) < 0.001)
		if(abs(lastError-error) < 10)
		{
			break;
		}
		lastError = error;

		/*for(long v=0; v<vectorCount; v++)
		{
			for(int c=0; c<k; c++)
			{
				forceMatrix->data[v][c] /= totalForCentre[c];
			}
		}*/

		/*for(int c=0; c<k; c++)
		{
			for(long f=0; f<featureCount; f++)
			{
				centres[c][f] = 0;
			}
			for(long v=0; v<vectorCount; v++)
			{
				for(long f=0; f<featureCount; f++)
				{
					//centres[c][f] += forceMatrix->data[v][c] * (*featureVectors)[v][f];
					centres[c][f] += (*forceMatrix)(v, c) * (*featureVectors)[v][f] / totalForCentre[c];
				}
			}
		}*/
		centres = (*forceMatrix).transpose() * (*featureVectors);
		for(int c=0; c<k; c++)
		{
			centres.row(c) /= forceMatrix->col(c).sum();
		}
	}

	return forceMatrix;
}