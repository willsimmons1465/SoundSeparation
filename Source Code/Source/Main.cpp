#include <iostream>
#include <string>
#include <vector>
#include "Cluster.h"
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

int main(void)
{
	//testFileManager();
	//testFFT();
	//testStft();
	//testComplex();
	//testHamming();
	//testCluster();
	testSinModSep();

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

	vector<vector<complex<double>>>* lstft = Transform::stft(&lSampleVector, 4096, 1024, 4096);
	vector<vector<complex<double>>>* rstft = Transform::stft(&rSampleVector, 4096, 1024, 4096);

	cout << "STFTs calculated" << endl;
	cout << lstft->size() << "x" << (*lstft)[0].size() << endl;
	cout << (*lstft)[0][0].real() << " + " << (*lstft)[0][0].imag() << "j" << endl;
	cout << (*lstft)[0][2048].real() << " + " << (*lstft)[0][2048].imag() << "j" << endl;

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
	WavFileManager fileManager("Guitar Trumpet/GD3vlfn_TGs41fn");

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
