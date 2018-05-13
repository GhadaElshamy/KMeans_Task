#include <iostream>
#include<omp.h>
#include<vector>
#include<fstream>
#include<string>
#include<sstream>
#include<stdlib.h>
#include<time.h>

using namespace std;

class datapoint
{
public:
	vector<float> coordinates;
	int cluster;
};



void AssignClusters(vector<datapoint>& data, vector<datapoint> clusters, vector<vector<datapoint>>& clusterPoints, int kcluster, int numPoints)
{
	int coord = data[0].coordinates.size();
	for (int i = 0; i < numPoints; i++)
	{
		float mindist = 10e9;
		int j;
		for (j = 0; j < kcluster; j++)
		{
			vector<float> cluster = clusters[j].coordinates;
			float dist = 0.0;
			for (int m = 0; m < coord; m++)
			{
				dist += (data[i].coordinates[m] - cluster[m])*(data[i].coordinates[m] - cluster[m]);
			}
			dist = sqrtf(dist);

			if (dist < mindist)
			{
				data[i].cluster = j;
				mindist = dist;
			}
		}
		clusterPoints[data[i].cluster].push_back(data[i]);
	}
}

vector<datapoint> CalculateNewClusters(vector<datapoint>& clusters, vector<vector<datapoint>> clusterPoints, int kcluster)
{
	vector<datapoint> newClusters;
	for (int i = 0; i < kcluster; i++)
	{
		int totalPoints = clusterPoints[i].size();
		int coord = clusters[0].coordinates.size();
		vector<float> mean;
		mean.assign(coord, 0);
		for (int j = 0; j < totalPoints; j++)
		{
			for (int m = 0; m < coord; m++)
			{
				mean[m] += (clusterPoints[i][j].coordinates[m])/totalPoints;
			}

		}
		datapoint obj = datapoint();
		obj.coordinates = mean;
		newClusters.push_back(obj);
	}

	return newClusters;
}

bool CheckClusters(vector<datapoint> oldClusters, vector<datapoint> newClusters,float threshold)
{
	for (int i = 0; i < oldClusters.size(); i++)
	{
		float dist = 0.0;
		for (int j = 0; j < oldClusters[i].coordinates.size(); j++)
		{
			dist += (oldClusters[i].coordinates[j] - newClusters[i].coordinates[j])*(oldClusters[i].coordinates[j] - newClusters[i].coordinates[j]);
		}
		dist = sqrtf(dist);

		if (dist > threshold)
			return false;
	}
	return true;
}

void ParallelAssignCluster(vector<datapoint>& data, vector<datapoint> clusters, vector<vector<datapoint>>& clusterPoints, int kcluster, int numPoints,int threadsNum)
{
	int coord = data[0].coordinates.size() ;
	float mindist = 10e9;
	float dist = 0.0;
	omp_set_num_threads(threadsNum); 
#pragma omp parallel
{
#pragma omp for firstprivate(dist) private(mindist)
	for (int i = 0; i < numPoints*kcluster*coord; i++) // same as using collapse(2) but built in collabse doesn't work right here
	{
		int n = i / (coord*kcluster), j = (i % (coord*kcluster)) / coord, m = (i % (coord*kcluster)) % coord;

		vector<float> cluster = clusters[j].coordinates;

		if (m == 0)
		{
			dist = 0.0;
		}
		dist += (data[n].coordinates[m] - cluster[m])*(data[n].coordinates[m] - cluster[m]);

		if (m == coord - 1)
		{
			dist = sqrtf(dist);

			if (dist < mindist)
			{
				data[n].cluster = j;
				mindist = dist;
			}
		}

		if (i % (kcluster*coord) == 0)
		{
			clusterPoints[data[n].cluster].push_back(data[n]);
			mindist = 10e9;
		}
	}
}
}

vector<datapoint> ParallelCalculateNewClusters(vector<datapoint>& clusters, vector<vector<datapoint>> clusterPoints, int kcluster,int threadNum)
{
	vector<datapoint> newClusters;
	omp_set_num_threads(threadNum);
	omp_set_nested(threadNum);

#pragma omp parallel
{
#pragma omp for
	for (int i = 0; i < kcluster; i++)
	{
		int totalPoints = clusterPoints[i].size();
		int coord = clusters[0].coordinates.size();
		vector<float> mean;
		mean.assign(coord, 0);
#pragma omp parallel for
		for (int j = 0; j < totalPoints*coord; j++)
		{
			int n = j / coord, m = j% coord;

           #pragma omp atomic
			//{
			  mean[m] += (clusterPoints[i][n].coordinates[m]) / totalPoints;
			//}

		}

		datapoint obj = datapoint();
		obj.coordinates = mean;
		newClusters.push_back(obj);
	}
}
	return newClusters;
}

bool ParallelCheckClusters(vector<datapoint> oldClusters, vector<datapoint> newClusters, float threshold,int threadsNum)
{
	omp_set_num_threads(threadsNum);
	int coord = oldClusters[0].coordinates.size();
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < oldClusters.size()*coord; i++)   // same as using collapse(2) but built in collabse doesn't work right here
		{
			int n = i / coord, j = i%coord;
			float dist = 0.0;

			if (j == 0)
			{
				dist = 0.0;
			}

			dist += (oldClusters[n].coordinates[j] - newClusters[n].coordinates[j])*(oldClusters[n].coordinates[j] - newClusters[n].coordinates[j]);

			if (j == oldClusters[0].coordinates.size() - 1)
			{
				dist = sqrtf(dist);

				if (dist > threshold)
					return false;
			}
		}
	}
	return true;
}

vector<datapoint> ReadData(string path, int& kclusters, int& datasize)
{
	ifstream file;
	vector<datapoint> data;
	try
	{
		file.open(path);
		string line = "",z="";
		
		getline(file, line, '\n');
		
		for (int i = 0; i < line.size(); i++)
		{
			if (line[i] == ' ' && datasize!=0)
			{
				datasize = stoi(z);
				z.clear();
			}
			else
			z += line[i];
		}
		kclusters = stoi(z);
		

		while (getline(file, line, '\n'))
		{
			vector<float> tmp;
			z.clear();
			for (int i = 0; i < line.size(); i++)
			{
				if (line[i] == ',')
				{
					tmp.push_back(stof(z));
					
					z.clear();
				}
				else
					z += line[i];
			}
			if (!z.empty())
			{
				tmp.push_back(stof(z));
			}
			datapoint obj = datapoint();
			obj.coordinates = tmp;
			data.push_back(obj);
		}

		file.close();
	}
	catch (exception e)
	{
		printf(e.what());
	}
	return data;
}

void WriteData(vector<datapoint> clusters, string infileName)
{
	string outfile = infileName + "_cluster_centers.txt";
	ofstream file(outfile);

	try
	{
		int coord = clusters[0].coordinates.size();
		for (int i = 0; i < clusters.size(); i++)
		{
			file << (i+1) << ',';
			for (int j = 0; j < coord; j++)
			{
				if (j == coord - 1)
					file << clusters[i].coordinates[j] << '\n';
				else
					file << clusters[i].coordinates[j] << ',';
			}
		}
		file.close();
	}
	catch (exception e)
	{
		cout << e.what();
	}
}

void SequentialKMeans(string filename,vector<datapoint> data, vector<datapoint> clusters, vector<vector<datapoint>> clusterPoints, int kcluster, int numPoints, float threshold = 0.001)
{
	while (true)
	{
		AssignClusters(data, clusters, clusterPoints, kcluster, numPoints);
		vector<datapoint> newClusters = CalculateNewClusters(clusters, clusterPoints, kcluster);
		if (CheckClusters(clusters, newClusters, threshold))
		{
			break;
		}
		else
			clusters = newClusters;
	}
	WriteData(clusters, filename);
}

void ParallelKMeans(string filename, vector<datapoint> data, vector<datapoint> clusters, vector<vector<datapoint>> clusterPoints, int kcluster, int numPoints,int threadsNum=4, float threshold = 0.001)
{
	// Collapse is not working in the right way  so i collapsed the for loops manualy

	while (true)
	{
		ParallelAssignCluster(data, clusters, clusterPoints, kcluster, numPoints,threadsNum);
		vector<datapoint> newClusters = ParallelCalculateNewClusters(clusters, clusterPoints, kcluster,threadsNum);
		if (ParallelCheckClusters(clusters, newClusters, threshold,threadsNum))
		{
			break;
		}
		else
			clusters = newClusters;
	}
	WriteData(clusters, filename);
}

void main()
{
 
#pragma region initialization
	int kclusters,numpoints;
	int threadsNum = 4;
	float threshold = 0.001;
	string filename = "IrisDataset";
	string path = "F:\\ghada\\.Trashes\\Fourth_Year[CS]\\DC\\Tasks\\Task_1 K-Means\\"+filename+".txt";
	vector<datapoint> data;
	vector<datapoint> clusters;
	vector<vector<datapoint>> clusterPoints;

	data = ReadData(path, kclusters, numpoints);

	// initial random clusters take 1st 5 clusters as initial clusters
	vector<int> tmp;
	for (int i = 0; i < kclusters; i++)
	{
			clusters.push_back(data[i]);
	}

	for (int i = 0; i < kclusters; i++)
	{
		vector<datapoint> tmp = {};
		clusterPoints.push_back(tmp);
	}
#pragma endregion


//SequentialKMeans(filename, data, clusters, clusterPoints, kclusters, numpoints, threshold);
ParallelKMeans(filename, data, clusters, clusterPoints, kclusters, numpoints, threadsNum, threshold);


}
