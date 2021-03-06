#include <iostream>
#include <string>
#include <vector>
#include "Cluster.h"
#include "NMFSeparation.h"
#include "SinusoidalModelSeparation.h"
#include "Transform.h"
#include "WavFileManager.h"

using namespace std;

void testFileManager(void);
void testFFT(void);
void testStft(void);
void testComplex(void);
void testHamming(void);
void testCluster(void);
void testSinModSep(void);
void testNMF(void);
void testNMFSep(void);
void testSoftCluster(void);

int main(void)
{
	//testFileManager();
	//testFFT();
	//testStft();
	//testComplex();
	//testHamming();
	//testCluster();
	//testSinModSep();
	//testNMF();
	//testNMFSep();
	//SinusoidalModelSeparation::optimiseParams();
	testSoftCluster();

	string wait;
	cin >> wait;

    return 0;
}

void testFileManager(void)
{
	WavFileManager fileManager("Guitar Trumpet/GD3vlfn_TGs41fn");

	vector<double> lSampleVector;
	vector<double> rSampleVector;
	fileManager.readSoundSample(lSampleVector, rSampleVector);

	cout << "Samples read" << endl;

	fileManager.writeDerivedOutput(0, lSampleVector, rSampleVector);

	cout << "Samples written" << endl;
}

void testFFT(void)
{
	double samples[] = {5, 2, 75, 86, -98, -64, 46, -26,
						94, 46, -41, -24, -32, -18, -54, 39,
						93, 76};
	vector<complex<double>>* hfft = Transform::halfFFT(samples, sizeof(samples)/sizeof(*samples));
	const int halfSize = hfft->size();
	for(int i=0; i<halfSize; i++)
	{
		cout << (*hfft)[i].real() << "+j" << (*hfft)[i].imag() << endl;
	}
	vector<complex<double>>* expandedFFT = new vector<complex<double>>();
	for(int i=0; i<halfSize; i++)
	{
		expandedFFT->push_back((*hfft)[i]);
	}
	const int eLength = 2*(hfft->size()-1);
	for(int i=halfSize; i<eLength; i++)
	{
		expandedFFT->push_back(conj((*hfft)[eLength-i]));
	}
	double* reconstructed = Transform::ifft(expandedFFT);
	for(int i=0; i<eLength; i++)
	{
		cout << reconstructed[i] << endl;
	}
	delete expandedFFT;
	delete reconstructed;
}

void testStft(void)
{
	WavFileManager fileManager("Guitar Trumpet/GD3vlfn_TGs41fn");

	vector<double> lSampleVector;
	vector<double> rSampleVector;
	fileManager.readSoundSample(lSampleVector, rSampleVector);

	cout << "Samples read" << endl;

	Matrix<complex<double>>* lstft = Transform::stft(&lSampleVector, 4096, 1024, 4096);
	Matrix<complex<double>>* rstft = Transform::stft(&rSampleVector, 4096, 1024, 4096);

	cout << "STFTs calculated" << endl;
	cout << lstft->data.size() << "x" << lstft->data[0].size() << endl;
	cout << lstft->data[0][0].real() << " + " << lstft->data[0][0].imag() << "j" << endl;
	cout << lstft->data[0][2048].real() << " + " << lstft->data[0][2048].imag() << "j" << endl;

	vector<double>* lReconstructed = Transform::istft(lstft, 1024);
	vector<double>* rReconstructed = Transform::istft(rstft, 1024);

	cout << "iSTFTs calculated" << endl;

	fileManager.writeDerivedOutput(0, *lReconstructed, *rReconstructed);

	cout << "Samples written" << endl;

	delete lstft;
	delete rstft;
	delete lReconstructed;
	delete rReconstructed;
}

void testComplex(void)
{
	complex<double> z = complex<double>(2., 1.) + complex<double>(2., 3.);
	cout << z.real() << " + " << z.imag() << "j" <<endl;
}

void testHamming(void)
{
	const long length = 40;
	double* window = Transform::hamming(length);
	for(long i=0; i<length; i++)
	{
		cout << window[i] << endl;
	}
	delete window;
}

void testCluster(void)
{
	srand(19873495);
	vector<vector<double>> featureVectors = vector<vector<double>>(10);
	featureVectors[0].push_back(0.);
	featureVectors[1].push_back(1.);
	featureVectors[2].push_back(2.);
	featureVectors[3].push_back(3.);
	featureVectors[4].push_back(4.);
	featureVectors[5].push_back(15.);
	featureVectors[6].push_back(16.);
	featureVectors[7].push_back(17.);
	featureVectors[8].push_back(18.);
	featureVectors[9].push_back(19.);
	for(int i=0; i<featureVectors.size(); i++)
	{
		featureVectors[i].push_back(16.-0.4*featureVectors[i][0]);
	}
	vector<int>* clusterTags = Cluster::kClusterLloyd(&featureVectors, 3);
	for(int i=0; i<clusterTags->size(); i++)
	{
		cout << (*clusterTags)[i] << endl;
	}
	delete clusterTags;
}

void testSinModSep(void)
{
	//WavFileManager fileManager("Guitar Trumpet/GD3vlfn_TGs41fn");
	WavFileManager fileManager("OptimiseTests/SE51fn_TGs41fn");

	vector<double> lSampleVector;
	vector<double> rSampleVector;
	fileManager.readSoundSample(lSampleVector, rSampleVector);

	cout << "Samples read" << endl;

	vector<vector<vector<double>>>* separated = SinusoidalModelSeparation::separate(&lSampleVector, &rSampleVector, 2);

	cout << "Separated" << endl;

	for(int i=0; i<2; i++)
	{
		fileManager.writeDerivedOutput(i, (*separated)[i][0], (*separated)[i][1]);
	}

	cout << "Separated sources written" << endl;

	delete separated;
}

void testNMF(void)
{
	Matrix<double> realMix0 = Matrix<double>(4,1);
	Matrix<double> realSource0 = Matrix<double>(1,4);
	realMix0.data[0][0] = 2.;
	realMix0.data[1][0] = 5.;
	realMix0.data[2][0] = 1.;
	realMix0.data[3][0] = 2.4;
	realSource0.data[0][0] = 3;
	realSource0.data[0][1] = 1;
	realSource0.data[0][2] = 4;
	realSource0.data[0][3] = 2;

	Matrix<double>* mat0 = Matrix<double>::multiply(&realMix0, &realSource0);

	Matrix<double> realMix1 = Matrix<double>(4,1);
	Matrix<double> realSource1 = Matrix<double>(1,4);
	realMix1.data[0][0] = 1.;
	realMix1.data[1][0] = 7.;
	realMix1.data[2][0] = 3.;
	realMix1.data[3][0] = 4.;
	realSource1.data[0][0] = 5;
	realSource1.data[0][1] = 1;
	realSource1.data[0][2] = 3;
	realSource1.data[0][3] = 2;

	Matrix<double>* mat1 = Matrix<double>::multiply(&realMix1, &realSource1);

	Matrix<double> mat = Matrix<double>(4,4);

	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			mat.data[i][j] = mat0->data[i][j] + mat1->data[i][j];
			cout << mat.data[i][j] << " ";
		}
		cout << endl;
	}

	cout << "Matrix created" << endl;

	vector<Matrix<double>>* nmfResult = Transform::nmf(&mat, 2);

	cout << "NMF completed" << endl;

	for(int i=0; i<4; i++)
	{
		cout << nmfResult->data()[0].data[i][0] << " ";
	}
	cout << endl;
	for(int i=0; i<4; i++)
	{
		cout << nmfResult->data()[1].data[0][i] << " ";
	}
	cout << endl;
	for(int i=0; i<4; i++)
	{
		cout << nmfResult->data()[0].data[i][1] << " ";
	}
	cout << endl;
	for(int i=0; i<4; i++)
	{
		cout << nmfResult->data()[1].data[1][i] << " ";
	}

	delete mat0;
	delete mat1;
	delete nmfResult;
}

void testNMFSep(void)
{
	//WavFileManager fileManager("Guitar Trumpet/GD3vlfn_TGs41fn");
	WavFileManager fileManager("OptimiseTests/SE51fn_TGs41fn");

	vector<double> lSampleVector;
	vector<double> rSampleVector;
	fileManager.readSoundSample(lSampleVector, rSampleVector);

	cout << "Samples read" << endl;

	vector<vector<vector<double>>>* separated = NMFSeparation::separate(&lSampleVector, &rSampleVector, 2, 1000);

	cout << "Separated" << endl;

	for(int i=0; i<2; i++)
	{
		fileManager.writeDerivedOutput(i, (*separated)[i][0], (*separated)[i][1]);
	}

	cout << "Separated sources written" << endl;

	delete separated;
}

void testSoftCluster(void)
{
	srand(19873495);
	vector<vector<double>> featureVectors = vector<vector<double>>(10);
	featureVectors[0].push_back(0.);
	featureVectors[1].push_back(1.);
	featureVectors[2].push_back(2.);
	featureVectors[3].push_back(3.);
	featureVectors[4].push_back(4.);
	featureVectors[5].push_back(15.);
	featureVectors[6].push_back(16.);
	featureVectors[7].push_back(17.);
	featureVectors[8].push_back(18.);
	featureVectors[9].push_back(19.);
	for(int i=0; i<featureVectors.size(); i++)
	{
		featureVectors[i].push_back(16.-0.4*featureVectors[i][0]);
	}
	Matrix<double>* clusters = Cluster::softKCluster(&featureVectors, 3, 1.);
	for(int i=0; i<clusters->data.size(); i++)
	{
		for(int j=0; j<clusters->data[0].size(); j++)
		{
			cout << clusters->data[i][j] << ", ";
		}
		cout << endl;
	}
	delete clusters;
}