#include <iostream>
#include <string>
#include <vector>
#include "BASSEval.h"
#include "Cluster.h"
//#include "NMFSeparation.h"
//#include "SinusoidalModelSeparation.h"
#include "Separation.h"
#include "SinusoidalTrajectory.h"
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

/*int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		cout << "Too few arguments" << endl;
		cout << "Usage: filename k <flags>" << endl;
		string wait;
		cin >> wait;
		return -1;
	}
	juce::String filename = argv[1];
	juce::String kstr = argv[2];
	int k = kstr.getIntValue();
	bool sinFlag = false;
	bool nmfFlag = false;
	bool hardFlag = false;
	bool softFlag = false;
	bool matFlag = false;
	bool naivFlag = false;
	bool revFlag = false;
	bool verbFlag = false;
	for (int i = 3; i < argc; i++)
	{
		juce::String flag = juce::String(argv[i]);
		if (flag.equalsIgnoreCase("-sin"))
		{
			sinFlag = true;
		}
		else if (flag.equalsIgnoreCase("-nmf"))
		{
			nmfFlag = true;
		}
		else if (flag.equalsIgnoreCase("-hc"))
		{
			hardFlag = true;
		}
		else if (flag.equalsIgnoreCase("-sc"))
		{
			softFlag = true;
		}
		else if (flag.equalsIgnoreCase("-mc"))
		{
			matFlag = true;
		}
		else if (flag.equalsIgnoreCase("-nc"))
		{
			naivFlag = true;
		}
		else if (flag.equalsIgnoreCase("-r"))
		{
			revFlag = true;
		}
		else if (flag.equalsIgnoreCase("-v"))
		{
			verbFlag = true;
		}
		else
		{
			cout << "Did not recognise: " << argv[i] << endl;
			cout << "Usage: filename k <flags>" << endl;
			string wait;
			cin >> wait;
			return -1;
		}
	}
	if (sinFlag && nmfFlag)
	{
		cout << "-sin and -nmf are mutually exclusive" << endl;
		string wait;
		cin >> wait;
		return -1;
	}
	if (hardFlag && (softFlag || matFlag || naivFlag)
		|| softFlag && (matFlag || naivFlag)
		|| matFlag && naivFlag)
	{
		cout << "-sc, -hc, -mc and -nc are mutually exclusive" << endl;
		string wait;
		cin >> wait;
		return -1;
	}
	if (nmfFlag && naivFlag)
	{
		cout << "Naive clustering only available with Sinusoidal Trajectories" << endl;
		string wait;
		cin >> wait;
		return -1;
	}

	Separation::FeatureOption feature = Separation::SINUSOIDS;
	if (nmfFlag)
	{
		feature = Separation::MATRIXFACTORS;
	}

	Separation::ClusterOption cluster = Separation::HARD;
	if (softFlag)
	{
		cluster = Separation::SOFT;
	}
	else if (matFlag)
	{
		cluster = Separation::MATRIX;
	}
	else if (naivFlag)
	{
		cluster = Separation::NAIVE;
	}

	WavFileManager fileManager(filename);

	srand(19873496);

	Eigen::VectorXd lSampleVector;
	Eigen::VectorXd rSampleVector;
	fileManager.readSoundSample(&lSampleVector, &rSampleVector);

	vector<vector<Eigen::VectorXd>>* separated = Separation::separate(&lSampleVector, &rSampleVector, k, feature, cluster, revFlag, verbFlag);

	for(int i=0; i<k; i++)
	{
		fileManager.writeDerivedOutput(i, &((*separated)[i][0]), &((*separated)[i][1]));
	}

	delete separated;

	string wait;
	cin >> wait;

    return 0;
}*/

/*int main(void)
{
	vector<WavFileManager*>* testSet = BASSEval::validationSetSounds();
	BASSEval::pairWiseStereoTest(testSet, Separation::MATRIXFACTORS, Separation::MATRIX, true, true);
	for(int i=0; i<testSet->size(); i++)
	{
		delete testSet->data()[i];
	}
	delete testSet;

	string wait;
	cin >> wait;

    return 0;
}*/

int main(void)
{
	BASSEval::standardTest(Separation::SINUSOIDS, Separation::SOFT, false);
	//BASSEval::incrementKTest(Separation::MATRIXFACTORS);
	//BASSEval::incrementStereoTest(Separation::SINUSOIDS);
	//BASSEval::incrementOffsetTest(Separation::SINUSOIDS);
	//BASSEval::incrementFrequencyTest(Separation::MATRIXFACTORS);
	//BASSEval::incrementPitchTest(Separation::MATRIXFACTORS);
	//BASSEval::noiseTest(Separation::SINUSOIDS, 2.05);

	string wait;
	cin >> wait;

	return 0;
}

/*int main(void)
{
	vector<WavFileManager*>* testSet = BASSEval::validationSetSounds();
	for(double val=50; val<=150; val += 10)
	{
		SinusoidalTrajectory::distanceMissPenaltyStereo = val;
		//Separation::softClusterStiffness = val;
		//Separation::matrixWeightMix = val;
		cout << "MPS = " << val << "; ";
		BASSEval::pairWiseStereoTest(testSet, Separation::SINUSOIDS, Separation::SOFT, false, false);
	}
	for(int i=0; i<testSet->size(); i++)
	{
		delete testSet->data()[i];
	}
	delete testSet;

	string wait;
	cin >> wait;

    return 0;
}*/

/*int main(void)
{
	vector<WavFileManager*>* testSet = BASSEval::testSetSounds();
	int i=4;
	int j=12;
	double stereoSep = 0.13;

	srand(19873496);

	Eigen::VectorXd li;
	Eigen::VectorXd ri;
	Eigen::VectorXd lj;
	Eigen::VectorXd rj;
	testSet->data()[i-1]->readSoundSample(&li, &ri);
	testSet->data()[j-1]->readSoundSample(&lj, &rj);

	double iStereo = 0.5 * (1 + stereoSep);
	double jStereo = 0.5 * (1 - stereoSep);

	long sampleCount = max(li.size(), lj.size());
	Eigen::VectorXd zeroVector = Eigen::VectorXd::Zero(sampleCount);

	Eigen::VectorXd iSampleVector = (li + ri)/sqrt(2.);
	Eigen::VectorXd jSampleVector = (lj + rj)/sqrt(2.);
	iSampleVector.conservativeResizeLike(zeroVector);
	jSampleVector.conservativeResizeLike(zeroVector);
	Eigen::VectorXd lSampleVector = sqrt(1-iStereo) * iSampleVector + sqrt(1-jStereo) * jSampleVector;
	Eigen::VectorXd rSampleVector = sqrt(iStereo) * iSampleVector + sqrt(jStereo) * jSampleVector;

	vector<vector<Eigen::VectorXd>>* separated = Separation::separate(&lSampleVector, &rSampleVector, 2, Separation::SINUSOIDS, Separation::SOFT, false, true);

	WavFileManager fileManager("C:\\Users\\Will\\OneDrive\\Uni\\Individual Project\\SoundSeparation\\OutputSpectrograms\\testOutput0.1.wav");

	for(int o=0; o<2; o++)
	{
		fileManager.writeDerivedOutput(i*1000 + j*10 + o, &((*separated)[o][0]), &((*separated)[o][1]));
	}


	for(int i=0; i<testSet->size(); i++)
	{
		delete testSet->data()[i];
	}
	delete testSet;
	delete separated;

	string wait;
	cin >> wait;

	return 0;
}*/

/*int main(void)
{
	//testFileManager();
	//testFFT();
	testStft();
	//testComplex();
	//testHamming();
	//testCluster();
	//testSinModSep();
	//testNMF();
	//testNMFSep();
	//SinusoidalModelSeparation::optimiseParams();
	//testSoftCluster();

	string wait;
	cin >> wait;

    return 0;
}*/

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
	double s[] = {5, 2, 75, 86, -98, -64, 46, -26,
						94, 46, -41, -24, -32, -18, -54, 39,
						93, 76};
	vector<double> samples = vector<double>(s, s + sizeof(s) / sizeof(*s));
	vector<complex<double>>* hfft = Transform::halfFFT(&samples);
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
	vector<double>* reconstructed = Transform::ifft(expandedFFT);
	for(int i=0; i<eLength; i++)
	{
		cout << (*reconstructed)[i] << endl;
	}
	delete hfft;
	delete expandedFFT;
	delete reconstructed;
}

void testStft(void)
{
	//WavFileManager fileManager("Guitar Trumpet/GD3vlfn_TGs41fn");
	WavFileManager fileManager("C:\\Users\\Will\\OneDrive\\Uni\\Individual Project\\SoundSeparation\\Sound Samples\\OptimiseTests\\SE51fn_TGs41fn.wav");
	
	vector<double> lSampleVector;
	vector<double> rSampleVector;
	fileManager.readSoundSample(lSampleVector, rSampleVector);

	cout << "Samples read" << endl;

	Eigen::MatrixXcd* lstft = Transform::stft(&lSampleVector, 4096, 1024, 4096);
	Eigen::MatrixXcd* rstft = Transform::stft(&rSampleVector, 4096, 1024, 4096);

	cout << "STFTs calculated" << endl;
	cout << lstft->rows() << "x" << lstft->cols() << endl;
	cout << (*lstft)(0, 0).real() << " + " << (*lstft)(0, 0).imag() << "j" << endl;
	cout << (*lstft)(0, 2048).real() << " + " << (*lstft)(0, 0).imag() << "j" << endl;

	Eigen::VectorXd* lVector = Transform::istft(lstft, 1024);
	Eigen::VectorXd* rVector = Transform::istft(rstft, 1024);
	vector<double> lReconstructed = vector<double>(lVector->data(), lVector->data() + lVector->size());
	vector<double> rReconstructed = vector<double>(rVector->data(), rVector->data() + rVector->size());

	cout << "iSTFTs calculated" << endl;

	fileManager.writeDerivedOutput(0, lReconstructed, rReconstructed);

	cout << "Samples written" << endl;

	delete lstft;
	delete rstft;
	delete lVector;
	delete rVector;
}

void testComplex(void)
{
	complex<double> z = complex<double>(2., 1.) + complex<double>(2., 3.);
	cout << z.real() << " + " << z.imag() << "j" <<endl;
}

void testHamming(void)
{
	const long length = 40;
	vector<double>* window = Transform::hamming(length);
	for(long i=0; i<length; i++)
	{
		cout << (*window)[i] << endl;
	}
	delete window;
}

void testCluster(void)
{
	srand(19873495);
	Eigen::MatrixXd featureVectors = Eigen::MatrixXd(10, 2);
	/*featureVectors[0].push_back(0.);
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
	}*/
	featureVectors(0, 0) = 0.;
	featureVectors(1, 0) = 1.;
	featureVectors(2, 0) = 2.;
	featureVectors(3, 0) = 3.;
	featureVectors(4, 0) = 4.;
	featureVectors(5, 0) = 15.;
	featureVectors(6, 0) = 16.;
	featureVectors(7, 0) = 17.;
	featureVectors(8, 0) = 18.;
	featureVectors(9, 0) = 19.;
	for(int i=0; i<featureVectors.rows(); i++)
	{
		featureVectors(i, 1) = 16.-0.4*featureVectors(i, 0);
	}
	vector<int>* clusterTags = Cluster::kClusterLloyd(&featureVectors, 3);
	for(int i=0; i<clusterTags->size(); i++)
	{
		cout << (*clusterTags)[i] << endl;
	}
	delete clusterTags;
}

/*void testSinModSep(void)
{
	//WavFileManager fileManager("Guitar Trumpet/GD3vlfn_TGs41fn");
	WavFileManager fileManager("OptimiseTests/SE51fn_TGs41fn");

	srand(19873496);

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
}*/

void testNMF(void)
{
	Eigen::MatrixXd realMix0 = Eigen::MatrixXd(4,1);
	Eigen::MatrixXd realSource0 = Eigen::MatrixXd(1,4);
	realMix0(0, 0) = 2.;
	realMix0(1, 0) = 5.;
	realMix0(2, 0) = 1.;
	realMix0(3, 0) = 2.4;
	realSource0(0, 0) = 3;
	realSource0(0, 1) = 1;
	realSource0(0, 2) = 4;
	realSource0(0, 3) = 2;

	Eigen::MatrixXd mat0 = realMix0 * realSource0;

	Eigen::MatrixXd realMix1 = Eigen::MatrixXd(4,1);
	Eigen::MatrixXd realSource1 = Eigen::MatrixXd(1,4);
	realMix1(0, 0) = 1.;
	realMix1(1, 0) = 7.;
	realMix1(2, 0) = 3.;
	realMix1(3, 0) = 4.;
	realSource1(0, 0) = 5;
	realSource1(0, 1) = 1;
	realSource1(0, 2) = 3;
	realSource1(0, 3) = 2;

	Eigen::MatrixXd mat1 = realMix1 * realSource1;

	//Eigen::MatrixXd mat = Eigen::MatrixXd(4,4);
	Eigen::MatrixXd mat = mat0 + mat1;

	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			//mat.data[i][j] = mat0->data[i][j] + mat1->data[i][j];
			cout << mat(i, j) << " ";
		}
		cout << endl;
	}

	cout << "Matrix created" << endl;

	vector<Eigen::MatrixXd>* nmfResult = Transform::nmf(&mat, 2);

	cout << "NMF completed" << endl;

	for(int i=0; i<4; i++)
	{
		cout << nmfResult->data()[0](i, 0) << " ";
	}
	cout << endl;
	for(int i=0; i<4; i++)
	{
		cout << nmfResult->data()[1](0, i) << " ";
	}
	cout << endl;
	for(int i=0; i<4; i++)
	{
		cout << nmfResult->data()[0](i, 1) << " ";
	}
	cout << endl;
	for(int i=0; i<4; i++)
	{
		cout << nmfResult->data()[1](1, i) << " ";
	}

	//delete mat0;
	//delete mat1;
	delete nmfResult;
}

/*void testNMFSep(void)
{
	//WavFileManager fileManager("Guitar Trumpet/GD3vlfn_TGs41fn");
	WavFileManager fileManager("OptimiseTests/SE51fn_TGs41fn");

	vector<double> lSampleVector;
	vector<double> rSampleVector;
	fileManager.readSoundSample(lSampleVector, rSampleVector);

	cout << "Samples read" << endl;

	vector<vector<vector<double>>>* separated = NMFSeparation::separate(&lSampleVector, &rSampleVector, 2, 10);

	cout << "Separated" << endl;

	for(int i=0; i<2; i++)
	{
		fileManager.writeDerivedOutput(i, (*separated)[i][0], (*separated)[i][1]);
	}

	cout << "Separated sources written" << endl;

	delete separated;
}*/

void testSoftCluster(void)
{
	srand(19873495);
	Eigen::MatrixXd featureVectors = Eigen::MatrixXd(10, 2);
	/*featureVectors[0].push_back(0.);
	featureVectors[1].push_back(1.);
	featureVectors[2].push_back(2.);
	featureVectors[3].push_back(3.);
	featureVectors[4].push_back(4.);
	featureVectors[5].push_back(15.);
	featureVectors[6].push_back(16.);
	featureVectors[7].push_back(17.);
	featureVectors[8].push_back(18.);
	featureVectors[9].push_back(19.);*/
	featureVectors(0, 0) = 0.;
	featureVectors(1, 0) = 1.;
	featureVectors(2, 0) = 2.;
	featureVectors(3, 0) = 3.;
	featureVectors(4, 0) = 4.;
	featureVectors(5, 0) = 15.;
	featureVectors(6, 0) = 16.;
	featureVectors(7, 0) = 17.;
	featureVectors(8, 0) = 18.;
	featureVectors(9, 0) = 19.;
	for(int i=0; i<featureVectors.rows(); i++)
	{
		featureVectors(i, 1) = 16.-0.4*featureVectors(i, 0);
	}
	Eigen::MatrixXd* clusters = Cluster::softKCluster(&featureVectors, 3, 1.);
	for(int i=0; i<clusters->rows(); i++)
	{
		for(int j=0; j<clusters->cols(); j++)
		{
			cout << (*clusters)(i, j) << ", ";
		}
		cout << endl;
	}
	delete clusters;
}