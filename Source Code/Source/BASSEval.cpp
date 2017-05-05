#include "BASSEval.h"
#include <algorithm>
#include <iostream>
#include <Eigen/LU>
#include <limits>

vector<WavFileManager*>* BASSEval::validationSetSounds()
{
	WavFileManager* guitar = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OptimiseTests/guitar_D3_very-long_forte_normal.wav");
	WavFileManager* trumpet = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OptimiseTests/trumpet_Gs4_1_forte_normal.wav");
	WavFileManager* clarinet = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OptimiseTests/clarinet_As3_1_forte_normal.wav");
	WavFileManager* saxophone = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OptimiseTests/saxophone_E5_1_forte_normal.wav");
	WavFileManager* violin = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OptimiseTests/violin_Fs6_1_forte_arco-normal.wav");
	vector<WavFileManager*>* sampleFiles = new vector<WavFileManager*>();
	sampleFiles->push_back(guitar);
	sampleFiles->push_back(trumpet);
	sampleFiles->push_back(clarinet);
	sampleFiles->push_back(saxophone);
	sampleFiles->push_back(violin);

	return sampleFiles;
}

vector<WavFileManager*>* BASSEval::testSetSounds()
{
	WavFileManager* clarinet0 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/clarinet_As5_05_pianissimo_normal.wav");
	WavFileManager* clarinet1 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/clarinet_B6_1_pianissimo_normal.wav");
	WavFileManager* clarinet2 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/clarinet_G4_1_piano_normal.wav");
	WavFileManager* clarinet3 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/clarinet_Gs5_025_piano_normal.wav");
	WavFileManager* guitar0 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/guitar_A4_very-long_piano_harmonics.wav");
	WavFileManager* guitar1 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/guitar_B5_very-long_piano_harmonics.wav");
	WavFileManager* guitar2 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/guitar_Cs4_very-long_forte_normal.wav");
	WavFileManager* guitar3 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/guitar_F2_very-long_piano_normal.wav");
	WavFileManager* saxophone0 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/saxophone_A4_15_forte_normal.wav");
	WavFileManager* saxophone1 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/saxophone_As3_025_pianissimo_normal.wav");
	WavFileManager* saxophone2 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/saxophone_D4_1_forte_normal.wav");
	WavFileManager* saxophone3 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/saxophone_Fs4_025_piano_normal.wav");
	WavFileManager* trumpet0 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/trumpet_C5_very-long_fortissimo_normal.wav");
	WavFileManager* trumpet1 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/trumpet_D5_05_fortissimo_normal.wav");
	WavFileManager* trumpet2 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/trumpet_E2_long_pianissimo_normal.wav");
	WavFileManager* trumpet3 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/trumpet_Fs4_long_piano_normal.wav");
	WavFileManager* violin0 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/violin_C4_025_forte_arco-normal.wav");
	WavFileManager* violin1 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/violin_Ds5_1_mezzo-forte_arco-normal.wav");
	WavFileManager* violin2 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/violin_E7_1_piano_arco-normal.wav");
	WavFileManager* violin3 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/violin_Gs7_1_mezzo-forte_natural-harmonic.wav");
	vector<WavFileManager*>* sampleFiles = new vector<WavFileManager*>();
	sampleFiles->push_back(clarinet0);
	sampleFiles->push_back(clarinet1);
	sampleFiles->push_back(clarinet2);
	sampleFiles->push_back(clarinet3);
	sampleFiles->push_back(guitar0);
	sampleFiles->push_back(guitar1);
	sampleFiles->push_back(guitar2);
	sampleFiles->push_back(guitar3);
	sampleFiles->push_back(saxophone0);
	sampleFiles->push_back(saxophone1);
	sampleFiles->push_back(saxophone2);
	sampleFiles->push_back(saxophone3);
	sampleFiles->push_back(trumpet0);
	sampleFiles->push_back(trumpet1);
	sampleFiles->push_back(trumpet2);
	sampleFiles->push_back(trumpet3);
	sampleFiles->push_back(violin0);
	sampleFiles->push_back(violin1);
	sampleFiles->push_back(violin2);
	sampleFiles->push_back(violin3);

	return sampleFiles;
}

/*Eigen::MatrixXd* projectionMatrix(vector<Eigen::VectorXd>* separated, vector<Eigen::VectorXd>* originals)
{
	int numOfSounds = originals->size();
	Eigen::MatrixXd* retVal = new Eigen::MatrixXd(numOfSounds, numOfSounds);
	for(int i=0; i<numOfSounds; i++)
	{
		for(int j=0; j<numOfSounds; j++)
		{
			(*retVal)(i, j) = separated->data()[i].dot(originals->data()[j]);
		}
	}
	return retVal;
}

Eigen::MatrixXd* GramMatrix(vector<Eigen::VectorXd>* originals)
{
	int numOfSounds = originals->size();
	Eigen::MatrixXd* retVal = new Eigen::MatrixXd(numOfSounds, numOfSounds);
	for(int i=0; i<numOfSounds; i++)
	{
		for(int j=i; j<numOfSounds; j++)
		{
			(*retVal)(i, j) = (*retVal)(j, i) = originals->data()[i].dot(originals->data()[j]);
		}
	}
	return retVal;
}*/



BASSEval::BASSStats BASSEval::stereoTest(vector<vector<Eigen::VectorXd>>* originals, Separation::FeatureOption fOp, Separation::ClusterOption cOp, bool reversible)
{
	int numOfSounds = originals->size();
	long maxLength = 0;
	for(int i=0; i<numOfSounds; i++)
	{
		long sampleCount = originals->data()[i][0].size();
		if (maxLength < sampleCount)
		{
			maxLength = sampleCount;
		}
	}

	//Mix sounds
	Eigen::VectorXd mixlSamples = Eigen::VectorXd::Zero(maxLength);
	Eigen::VectorXd mixrSamples = Eigen::VectorXd::Zero(maxLength);
	Eigen::VectorXd zeroVector = Eigen::VectorXd::Zero(maxLength);
	for(int i=0; i<numOfSounds; i++)
	{
		originals->data()[i][0].conservativeResizeLike(zeroVector);
		originals->data()[i][1].conservativeResizeLike(zeroVector);
		mixlSamples += originals->data()[i][0];
		mixrSamples += originals->data()[i][1];
	}

	//Perform separation
	srand(19873495);
	vector<vector<Eigen::VectorXd>>* separated = Separation::separate(&mixlSamples, &mixrSamples, numOfSounds, fOp, cOp, reversible, false);
	long recLength = separated->data()[0][0].size();

	//Concatenate left and right channels
	Eigen::MatrixXd concatOriginals = Eigen::MatrixXd(2*recLength, numOfSounds);
	zeroVector = Eigen::VectorXd::Zero(recLength);
	for(int i=0; i<numOfSounds; i++)
	{
		Eigen::VectorXd concat = Eigen::VectorXd(2*recLength);
		originals->data()[i][0].conservativeResizeLike(zeroVector);
		originals->data()[i][1].conservativeResizeLike(zeroVector);
		concat << originals->data()[i][0], originals->data()[i][1];
		concatOriginals.col(i) = concat;
	}
	Eigen::MatrixXd concatSeparated = Eigen::MatrixXd(2*recLength, numOfSounds);
	for(int i=0; i<numOfSounds; i++)
	{
		Eigen::VectorXd concat = Eigen::VectorXd(2*recLength);
		concat << separated->data()[i][0], separated->data()[i][1];
		concatSeparated.col(i) = concat;
	}

	//concatSeparated = concatOriginals;

	//Gather statistics
	Eigen::MatrixXd Rss = concatOriginals.transpose() * concatOriginals;
	Eigen::MatrixXd projections = concatSeparated.transpose() * concatOriginals;
	Eigen::MatrixXd expansions = projections * Rss.inverse().conjugate().transpose() * concatOriginals.transpose(); //Row i gives expansion of separated i
	/*for(int i=0; i<numOfSounds; i++)
	{
		Eigen::VectorXd proj = concatSeparated.col(i).transpose() * concatOriginals;
		Eigen::VectorXd c = Rss.inverse() * proj;
		expansions.row(i).setZero();
		for(int j=0; j<numOfSounds; j++)
		{
			expansions.row(i) += c(j) * concatOriginals.col(j).transpose();
		}
	}*/
	Eigen::MatrixXd SDR = Eigen::MatrixXd(numOfSounds, numOfSounds);
	Eigen::MatrixXd SIR = Eigen::MatrixXd(numOfSounds, numOfSounds);
	Eigen::MatrixXd SAR = Eigen::MatrixXd(numOfSounds, numOfSounds);
	/*cout << "Rss:" << endl << Rss << endl;
	cout << "projections:" << endl << projections << endl;*/
	for(int i=0; i<numOfSounds; i++)
	{
		for(int j=0; j<numOfSounds; j++)
		{
			//WavFileManager* wav = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OptimiseTests/SE51fn_TGs41fn.wav");
			Eigen::VectorXd s_target = (projections(i, j)/Rss(j, j)) * concatOriginals.col(j);
			Eigen::VectorXd e_interf = expansions.row(i).transpose() - s_target;
			Eigen::VectorXd e_artef = concatSeparated.col(i) - expansions.row(i).transpose();
			/*Eigen::VectorXd s_targetl = s_target.head(recLength);
			Eigen::VectorXd s_targetr = s_target.tail(recLength);
			Eigen::VectorXd e_interfl = e_interf.head(recLength);
			Eigen::VectorXd e_interfr = e_interf.tail(recLength);
			Eigen::VectorXd e_artefl = e_artef.head(recLength);
			Eigen::VectorXd e_artefr = e_artef.tail(recLength);
			wav->writeDerivedOutput(4*i, &(separated->data()[i][0]), &(separated->data()[i][1]));
			//wav->writeDerivedOutput(4*i, &mixlSamples, &mixrSamples);
			wav->writeDerivedOutput(4*i + 1, &s_targetl, &s_targetr);
			wav->writeDerivedOutput(4*i + 2, &e_interfl, &e_interfr);
			wav->writeDerivedOutput(4*i + 3, &e_artefl, &e_artefr);
			cout << "proj/Rss == " << projections(i, j)/Rss(j, j) << endl;
			cout << "original.squaredNorm() == " << concatOriginals.col(j).squaredNorm() << endl;
			cout << "separated.squaredNorm() == " << concatSeparated.col(i).squaredNorm() << endl;
			cout << "s_target.squaredNorm() == " << s_target.squaredNorm() << endl;
			cout << "e_interf.squaredNorm() == " << e_interf.squaredNorm() << endl;
			cout << "e_artef.squaredNorm() == " << e_artef.squaredNorm() << endl;*/
			SDR(i, j) = 10 * log10(s_target.squaredNorm() / (e_interf + e_artef).squaredNorm());
			SIR(i, j) = 10 * log10(s_target.squaredNorm() / e_interf.squaredNorm());
			SAR(i, j) = 10 * log10((s_target + e_interf).squaredNorm() / e_artef.squaredNorm());

			//delete wav;
		}
	}

	double sumSDR = numeric_limits<double>::lowest();
	double sumSIR = 0;
	double sumSAR = 0;
	/*bool* mapped = new bool[numOfSounds];
	for(int i=0; i<numOfSounds; i++)
	{
		mapped[i] = false;
	}
	for(int i=0; i<numOfSounds; i++)
	{
		double maxSDR = numeric_limits<double>::lowest();
		int mapTo = -1;
		for(int j=0; j<numOfSounds; j++)
		{
			if(mapped[j])
			{
				continue;
			}
			double sdr = SDR(i, j);
			
			if(maxSDR < sdr)
			{
				maxSDR = sdr;
				mapTo = j;
			}
		}
		mapped[mapTo] = true;

		cout << mapTo;

		sumSDR += SDR(i, mapTo);
		sumSIR += SIR(i, mapTo);
		sumSAR += SAR(i, mapTo);
	}

	delete mapped;*/
	int* map = new int[numOfSounds];
	for(int i=0; i<numOfSounds; i++)
	{
		map[i] = i;
	}
	do {
		double sdr = 0.;
		double sir = 0.;
		double sar = 0.;
		for(int i=0; i<numOfSounds; i++)
		{
			int mapTo = map[i];
			sdr += SDR(i, mapTo);
			sir += SIR(i, mapTo);
			sar += SAR(i, mapTo);
		}
		if(sumSDR < sdr)
		{
			sumSDR = sdr;
			sumSIR = sir;
			sumSAR = sar;
		}
	} while (next_permutation(map, map+numOfSounds));
	delete map;

	/*cout << "SDR:" << endl << SDR << endl;
	cout << "SIR:" << endl << SIR << endl;
	cout << "SAR:" << endl << SAR << endl;*/
	/*cout << "SDR: " << sumSDR/numOfSounds << endl;
	cout << "SIR: " << sumSIR/numOfSounds << endl;
	cout << "SAR: " << sumSAR/numOfSounds << endl;*/

	delete separated;

	BASSStats retVal = BASSStats();
	retVal.avgSDR = sumSDR/numOfSounds;
	retVal.avgSIR = sumSIR/numOfSounds;
	retVal.avgSAR = sumSAR/numOfSounds;

	return retVal;
}

void BASSEval::pairWiseStereoTest(vector<WavFileManager*>* sampleFiles, Separation::FeatureOption fOp, Separation::ClusterOption cOp, bool reversible, bool fullResults)
{
	int numOfSounds = sampleFiles->size();
	vector<double> avgSDRs = vector<double>();
	vector<double> avgSIRs = vector<double>();
	vector<double> avgSARs = vector<double>();
	for(int i=0; i<numOfSounds; i++)
	{
		for(int j=i+1; j<numOfSounds; j++)
		{
			//Get samples from each sound
			Eigen::VectorXd ilSamples;
			Eigen::VectorXd irSamples;
			Eigen::VectorXd jlSamples;
			Eigen::VectorXd jrSamples;
			sampleFiles->data()[i]->readSoundSample(&ilSamples, &irSamples);
			sampleFiles->data()[j]->readSoundSample(&jlSamples, &jrSamples);

			//Generate stereo positions
			double stereoSplit = ((i+j) % numOfSounds) / (numOfSounds-1);
			double iStereo = 0.5 * (1 + stereoSplit);
			double jStereo = 0.5 * (1 - stereoSplit);

			//Pan sounds
			Eigen::VectorXd iSamples = (ilSamples + irSamples)/sqrt(2.);
			ilSamples = sqrt(1-iStereo) * iSamples;
			irSamples = sqrt(iStereo) * iSamples;
			Eigen::VectorXd jSamples = (jlSamples + jrSamples)/sqrt(2.);
			jlSamples = sqrt(1-jStereo) * jSamples;
			jrSamples = sqrt(jStereo) * jSamples;

			//Run test
			vector<vector<Eigen::VectorXd>> originalsToMix = vector<vector<Eigen::VectorXd>>(2);
			originalsToMix[0].push_back(ilSamples);
			originalsToMix[0].push_back(irSamples);
			originalsToMix[1].push_back(jlSamples);
			originalsToMix[1].push_back(jrSamples);
			//cout << i << ", " << j << endl;
			cout << "*";
			BASSStats results = stereoTest(&originalsToMix, fOp, cOp, reversible);
			avgSDRs.push_back(results.avgSDR);
			avgSIRs.push_back(results.avgSIR);
			avgSARs.push_back(results.avgSAR);
		}
	}

	cout << "SDR:" << endl;
	double sumSDR = 0;
	for(int i=0; i<avgSDRs.size(); i++)
	{
		if (fullResults)
		{
			cout << avgSDRs[i] << endl;
		}
		sumSDR += avgSDRs[i];
	}
	cout << "avg: " << sumSDR/avgSDRs.size() << endl;

	cout << "SIR:" << endl;
	double sumSIR = 0;
	for(int i=0; i<avgSIRs.size(); i++)
	{
		if (fullResults)
		{
			cout << avgSIRs[i] << endl;
		}
		sumSIR += avgSIRs[i];
	}
	cout << "avg: " << sumSIR/avgSIRs.size() << endl;

	cout << "SAR:" << endl;
	double sumSAR = 0;
	for(int i=0; i<avgSARs.size(); i++)
	{
		if (fullResults)
		{
			cout << avgSARs[i] << endl;
		}
		sumSAR += avgSARs[i];
	}
	cout << "avg: " << sumSAR/avgSARs.size() << endl;
}

void BASSEval::standardTest(Separation::FeatureOption fOp, Separation::ClusterOption cOp, bool reversible)
{
	vector<WavFileManager*>* sampleFiles = testSetSounds();

	/*int numOfTests = 10;
	int samples0[] = {20, 15, 14, 7, 1, 19, 3, 12, 17, 6};
	int samples1[] = {4, 8, 11, 16, 10, 5, 18, 13, 2, 9};
	double stereoSep[] = {0.02, 0.26, 0.73, 0.15, 0.42, 0.15, 0.83, 0.98, 0.53, 0.99};*/
	int numOfTests = 40;
	int samples0[] = {1, 1, 1, 1,
		2, 2, 2, 2, 
		3, 3, 3, 3, 
		4, 4, 4,
		5, 5, 5,
		6, 6, 6,
		7, 7, 7,
		8, 8, 8, 8,
		9, 9,
		10, 10,
		11,
		12,
		13, 13,
		14, 14, 14,
		16
	};
	int samples1[] = {11, 14, 17, 20,
		6, 11, 18, 19,
		4, 7, 9, 13,
		5, 9, 12,
		13, 16, 17,
		10, 12, 18,
		12, 15, 18,
		10, 15, 17, 18,
		19, 20,
		11, 17,
		15,
		20,
		16, 20,
		15, 16, 19,
		19
	};
	double stereoSep[] = {0.2, 0.29, 0.67, 0.53,
		0.6, 0.8, 0.17, 0.65,
		0.69, 0.95, 1., 0.8,
		0.59, 0.87, 0.13, 0.34,
		0.88, 0.33, 0.15, 0.02,
		0.69, 0.48, 0.75, 0.03,
		0.75, 0.74, 0.33, 0.86,
		0.47, 0.15, 0.23, 0.46,
		0.39, 0.4, 0.54, 0.17,
		0.66, 0.82, 0.8, 0.51
	};

	Eigen::VectorXd avgSDRs = Eigen::VectorXd(numOfTests);
	Eigen::VectorXd avgSIRs = Eigen::VectorXd(numOfTests);
	Eigen::VectorXd avgSARs = Eigen::VectorXd(numOfTests);
	for(int i=0; i<numOfTests; i++)
	{
		//Get samples from each sound
		Eigen::VectorXd ilSamples;
		Eigen::VectorXd irSamples;
		Eigen::VectorXd jlSamples;
		Eigen::VectorXd jrSamples;
		sampleFiles->data()[samples0[i]-1]->readSoundSample(&ilSamples, &irSamples);
		sampleFiles->data()[samples1[i]-1]->readSoundSample(&jlSamples, &jrSamples);
		
		//Generate stereo positions
		double iStereo = 0.5 * (1 + stereoSep[i]);
		double jStereo = 0.5 * (1 - stereoSep[i]);
		//double iStereo = 1.;
		//double jStereo = 0.;

		//Pan sounds
		Eigen::VectorXd iSamples = (ilSamples + irSamples)/sqrt(2.);
		ilSamples = sqrt(1-iStereo) * iSamples;
		irSamples = sqrt(iStereo) * iSamples;
		Eigen::VectorXd jSamples = (jlSamples + jrSamples)/sqrt(2.);
		jlSamples = sqrt(1-jStereo) * jSamples;
		jrSamples = sqrt(jStereo) * jSamples;

		//Run test
		vector<vector<Eigen::VectorXd>> originalsToMix = vector<vector<Eigen::VectorXd>>(2);
		originalsToMix[0].push_back(ilSamples);
		originalsToMix[0].push_back(irSamples);
		originalsToMix[1].push_back(jlSamples);
		originalsToMix[1].push_back(jrSamples);
		cout << "*";
		BASSStats results = stereoTest(&originalsToMix, fOp, cOp, reversible);
		avgSDRs(i) = results.avgSDR;
		avgSIRs(i) = results.avgSIR;
		avgSARs(i) = results.avgSAR;

		/*cout << i << endl;
		cout << "SDR " << results.avgSDR << endl;
		cout << "SIR " << results.avgSIR << endl;
		cout << "SAR " << results.avgSAR << endl;*/
	}

	cout << endl;
	double meanSDR = avgSDRs.sum() / numOfTests;
	double sigmaSDR = sqrt((avgSDRs.array() - meanSDR).square().sum() / (numOfTests-1));
	cout << "SDR: " << meanSDR << " +- " << sigmaSDR << endl;
	sort(avgSDRs.data(), avgSDRs.data() + avgSDRs.size());
	cout << avgSDRs(0) << ", " << avgSDRs(numOfTests/4) << ", " << avgSDRs(numOfTests/2) << ", " << avgSDRs(3*numOfTests/4) << ", " << avgSDRs(numOfTests-1) << endl;
	double meanSIR = avgSIRs.sum() / numOfTests;
	double sigmaSIR = sqrt((avgSIRs.array() - meanSIR).square().sum() / (numOfTests-1));
	cout << "SIR: " << meanSIR << " +- " << sigmaSIR << endl;
	sort(avgSIRs.data(), avgSIRs.data() + avgSIRs.size());
	cout << avgSIRs(0) << ", " << avgSIRs(numOfTests/4) << ", " << avgSIRs(numOfTests/2) << ", " << avgSIRs(3*numOfTests/4) << ", " << avgSIRs(numOfTests-1) << endl;
	double meanSAR = avgSARs.sum() / numOfTests;
	double sigmaSAR = sqrt((avgSARs.array() - meanSAR).square().sum() / (numOfTests-1));
	cout << "SAR: " << meanSAR << " +- " << sigmaSAR << endl;
	sort(avgSARs.data(), avgSARs.data() + avgSARs.size());
	cout << avgSARs(0) << ", " << avgSARs(numOfTests/4) << ", " << avgSARs(numOfTests/2) << ", " << avgSARs(3*numOfTests/4) << ", " << avgSARs(numOfTests-1) << endl;

	for(int i=0; i<sampleFiles->size(); i++)
	{
		delete sampleFiles->data()[i];
	}
	delete sampleFiles;
}

void BASSEval::incrementKTest(Separation::FeatureOption fOp)
{
	vector<WavFileManager*>* sampleFiles = testSetSounds();

	int samplesInTest[] = {2, 7, 8, 13, 18};
	//int samplesInTest[] = {2, 5, 10, 15, 18};
	double stereoPos[] = {0.93, 0.49, 0.39, 0.55, 0.32};
	//double stereoPos[] = {0.47, 0.15, 0.07, 0.68, 0.36};
	//double stereoPos[] = {0., 0.25, 0.5, 0.75, 1.};

	vector<vector<Eigen::VectorXd>> originalsToMix = vector<vector<Eigen::VectorXd>>(1);
	Eigen::VectorXd ilSamples;
	Eigen::VectorXd irSamples;
	sampleFiles->data()[samplesInTest[0]]->readSoundSample(&ilSamples, &irSamples);
	Eigen::VectorXd iSamples = (ilSamples + irSamples)/sqrt(2.);
	ilSamples = sqrt(1-stereoPos[0]) * iSamples;
	irSamples = sqrt(stereoPos[0]) * iSamples;
	originalsToMix[0].push_back(ilSamples);
	originalsToMix[0].push_back(irSamples);
	for(int k=2; k<=5; k++)
	{
		Eigen::VectorXd jlSamples;
		Eigen::VectorXd jrSamples;
		sampleFiles->data()[samplesInTest[k-1]]->readSoundSample(&jlSamples, &jrSamples);
		Eigen::VectorXd jSamples = (jlSamples + jrSamples)/sqrt(2.);
		jlSamples = sqrt(1-stereoPos[k-1]) * jSamples;
		jrSamples = sqrt(stereoPos[k-1]) * jSamples;
		originalsToMix.push_back(vector<Eigen::VectorXd>());
		originalsToMix[k-1].push_back(jlSamples);
		originalsToMix[k-1].push_back(jrSamples);

		BASSStats results = stereoTest(&originalsToMix, fOp, Separation::SOFT, false);

		cout << "k = " << k << endl;
		cout << "SDR = " << results.avgSDR << endl;
		cout << "SIR = " << results.avgSIR << endl;
		cout << "SAR = " << results.avgSAR << endl;
	}

	for(int i=0; i<sampleFiles->size(); i++)
	{
		delete sampleFiles->data()[i];
	}
	delete sampleFiles;
}

void BASSEval::incrementStereoTest(Separation::FeatureOption fOp)
{
	vector<WavFileManager*>* sampleFiles = testSetSounds();

	int samplesInTest[] = {7, 16};

	Eigen::VectorXd ilSamples;
	Eigen::VectorXd irSamples;
	Eigen::VectorXd jlSamples;
	Eigen::VectorXd jrSamples;
	sampleFiles->data()[samplesInTest[0]]->readSoundSample(&ilSamples, &irSamples);
	sampleFiles->data()[samplesInTest[1]]->readSoundSample(&jlSamples, &jrSamples);
	Eigen::VectorXd iSamples = (ilSamples + irSamples)/sqrt(2.);
	Eigen::VectorXd jSamples = (jlSamples + jrSamples)/sqrt(2.);
	for(double stereoPos=0.5; stereoPos<=1; stereoPos+=0.1)
	{
		ilSamples = sqrt(1-stereoPos) * iSamples;
		irSamples = sqrt(stereoPos) * iSamples;
		jlSamples = sqrt(stereoPos) * jSamples;
		jrSamples = sqrt(1-stereoPos) * jSamples;
		vector<vector<Eigen::VectorXd>> originalsToMix = vector<vector<Eigen::VectorXd>>(2);
		originalsToMix[0].push_back(ilSamples);
		originalsToMix[0].push_back(irSamples);
		originalsToMix[1].push_back(jlSamples);
		originalsToMix[1].push_back(jrSamples);

		BASSStats results = stereoTest(&originalsToMix, fOp, Separation::SOFT, false);

		cout << "stereoPos = " << stereoPos << endl;
		cout << "SDR = " << results.avgSDR << endl;
		cout << "SIR = " << results.avgSIR << endl;
		cout << "SAR = " << results.avgSAR << endl;
	}

	for(int i=0; i<sampleFiles->size(); i++)
	{
		delete sampleFiles->data()[i];
	}
	delete sampleFiles;
}

void BASSEval::incrementOffsetTest(Separation::FeatureOption fOp)
{
	WavFileManager* guitar = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OffsetTest/guitar_E4_very-long_forte_normal.wav");
	WavFileManager* s0 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OffsetTest/szero.wav");
	WavFileManager* s50 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OffsetTest/s50.wav");
	WavFileManager* s100 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OffsetTest/s100.wav");
	WavFileManager* s150 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OffsetTest/s150.wav");
	WavFileManager* s200 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/OffsetTest/s200.wav");
	vector<WavFileManager*>* sampleFiles = new vector<WavFileManager*>();
	sampleFiles->push_back(guitar);
	sampleFiles->push_back(s0);
	sampleFiles->push_back(s50);
	sampleFiles->push_back(s100);
	sampleFiles->push_back(s150);
	sampleFiles->push_back(s200);
	
	Eigen::VectorXd ilSamples;
	Eigen::VectorXd irSamples;
	sampleFiles->data()[0]->readSoundSample(&ilSamples, &irSamples);
	for(int i=0; i<5; i++)
	{
		Eigen::VectorXd jlSamples;
		Eigen::VectorXd jrSamples;
		sampleFiles->data()[i+1]->readSoundSample(&jlSamples, &jrSamples);
		vector<vector<Eigen::VectorXd>> originalsToMix = vector<vector<Eigen::VectorXd>>(2);
		originalsToMix[0].push_back(ilSamples);
		originalsToMix[0].push_back(irSamples);
		originalsToMix[1].push_back(jlSamples);
		originalsToMix[1].push_back(jrSamples);

		BASSStats results = stereoTest(&originalsToMix, fOp, Separation::SOFT, false);

		cout << "offset = " << 50*i << "ms" << endl;
		cout << "SDR = " << results.avgSDR << endl;
		cout << "SIR = " << results.avgSIR << endl;
		cout << "SAR = " << results.avgSAR << endl;
	}

	for(int i=0; i<sampleFiles->size(); i++)
	{
		delete sampleFiles->data()[i];
	}
	delete sampleFiles;
}

void BASSEval::incrementFrequencyTest(Separation::FeatureOption fOp)
{
	WavFileManager* base = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square465.wav");
	WavFileManager* s0 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/sin465.wav");
	WavFileManager* s10 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/sin475.wav");
	WavFileManager* s20 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/sin485.wav");
	WavFileManager* s30 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/sin495.wav");
	WavFileManager* s40 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/sin505.wav");
	WavFileManager* s50 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/sin515.wav");
	vector<WavFileManager*>* sampleFiles = new vector<WavFileManager*>();
	sampleFiles->push_back(base);
	sampleFiles->push_back(s0);
	sampleFiles->push_back(s10);
	sampleFiles->push_back(s20);
	sampleFiles->push_back(s30);
	sampleFiles->push_back(s40);
	sampleFiles->push_back(s50);
	
	Eigen::VectorXd ilSamples;
	Eigen::VectorXd irSamples;
	sampleFiles->data()[0]->readSoundSample(&ilSamples, &irSamples);
	for(int i=0; i<6; i++)
	{
		Eigen::VectorXd jlSamples;
		Eigen::VectorXd jrSamples;
		sampleFiles->data()[i+1]->readSoundSample(&jlSamples, &jrSamples);
		vector<vector<Eigen::VectorXd>> originalsToMix = vector<vector<Eigen::VectorXd>>(2);
		originalsToMix[0].push_back(ilSamples);
		originalsToMix[0].push_back(irSamples);
		originalsToMix[1].push_back(jlSamples);
		originalsToMix[1].push_back(jrSamples);

		BASSStats results = stereoTest(&originalsToMix, fOp, Separation::SOFT, false);

		cout << "offset = " << 10*i << "Hz" << endl;
		cout << "SDR = " << results.avgSDR << endl;
		cout << "SIR = " << results.avgSIR << endl;
		cout << "SAR = " << results.avgSAR << endl;
	}

	for(int i=0; i<sampleFiles->size(); i++)
	{
		delete sampleFiles->data()[i];
	}
	delete sampleFiles;
}

void BASSEval::incrementPitchTest(Separation::FeatureOption fOp)
{
	WavFileManager* base = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square465.wav");
	WavFileManager* s1 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square494.wav");
	WavFileManager* s2 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square523.wav");
	WavFileManager* s3 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square554.wav");
	WavFileManager* s4 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square587.wav");
	WavFileManager* s5 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square622.wav");
	WavFileManager* s6 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square659.wav");
	WavFileManager* s7 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square698.wav");
	WavFileManager* s8 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square740.wav");
	WavFileManager* s9 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square784.wav");
	WavFileManager* s10 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square831.wav");
	WavFileManager* s11 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square880.wav");
	WavFileManager* s12 = new WavFileManager("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/FreqTest/square932.wav");
	vector<WavFileManager*>* sampleFiles = new vector<WavFileManager*>();
	sampleFiles->push_back(base);
	sampleFiles->push_back(s1);
	sampleFiles->push_back(s2);
	sampleFiles->push_back(s3);
	sampleFiles->push_back(s4);
	sampleFiles->push_back(s5);
	sampleFiles->push_back(s6);
	sampleFiles->push_back(s7);
	sampleFiles->push_back(s8);
	sampleFiles->push_back(s9);
	sampleFiles->push_back(s10);
	sampleFiles->push_back(s11);
	sampleFiles->push_back(s12);
	
	Eigen::VectorXd ilSamples;
	Eigen::VectorXd irSamples;
	sampleFiles->data()[0]->readSoundSample(&ilSamples, &irSamples);
	for(int i=1; i<=12; i++)
	{
		Eigen::VectorXd jlSamples;
		Eigen::VectorXd jrSamples;
		sampleFiles->data()[i]->readSoundSample(&jlSamples, &jrSamples);
		vector<vector<Eigen::VectorXd>> originalsToMix = vector<vector<Eigen::VectorXd>>(2);
		originalsToMix[0].push_back(ilSamples);
		originalsToMix[0].push_back(irSamples);
		originalsToMix[1].push_back(jlSamples);
		originalsToMix[1].push_back(jrSamples);

		BASSStats results = stereoTest(&originalsToMix, fOp, Separation::SOFT, false);

		cout << "offset = " << i << endl;
		cout << "SDR = " << results.avgSDR << endl;
		cout << "SIR = " << results.avgSIR << endl;
		cout << "SAR = " << results.avgSAR << endl;
	}

	for(int i=0; i<sampleFiles->size(); i++)
	{
		delete sampleFiles->data()[i];
	}
	delete sampleFiles;
}

void BASSEval::noiseTest(Separation::FeatureOption fOp, double noiseMult)
{
	vector<WavFileManager*>* sampleFiles = testSetSounds();
	WavFileManager noiseFile("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/StandardTest/noise.wav");

	int numOfTests = 40;
	int samples0[] = {1, 1, 1, 1,
		2, 2, 2, 2, 
		3, 3, 3, 3, 
		4, 4, 4,
		5, 5, 5,
		6, 6, 6,
		7, 7, 7,
		8, 8, 8, 8,
		9, 9,
		10, 10,
		11,
		12,
		13, 13,
		14, 14, 14,
		16
	};
	int samples1[] = {11, 14, 17, 20,
		6, 11, 18, 19,
		4, 7, 9, 13,
		5, 9, 12,
		13, 16, 17,
		10, 12, 18,
		12, 15, 18,
		10, 15, 17, 18,
		19, 20,
		11, 17,
		15,
		20,
		16, 20,
		15, 16, 19,
		19
	};
	double stereoSep[] = {0.2, 0.29, 0.67, 0.53,
		0.6, 0.8, 0.17, 0.65,
		0.69, 0.95, 1., 0.8,
		0.59, 0.87, 0.13, 0.34,
		0.88, 0.33, 0.15, 0.02,
		0.69, 0.48, 0.75, 0.03,
		0.75, 0.74, 0.33, 0.86,
		0.47, 0.15, 0.23, 0.46,
		0.39, 0.4, 0.54, 0.17,
		0.66, 0.82, 0.8, 0.51
	};

	Eigen::VectorXd noiselSamples;
	Eigen::VectorXd noiserSamples;
	noiseFile.readSoundSample(&noiselSamples, &noiserSamples);
	noiselSamples = noiselSamples * noiseMult;
	noiserSamples = noiserSamples * noiseMult;

	Eigen::VectorXd inSNRs = Eigen::VectorXd(numOfTests);
	Eigen::VectorXd avgSDRs = Eigen::VectorXd(numOfTests);
	Eigen::VectorXd avgSIRs = Eigen::VectorXd(numOfTests);
	Eigen::VectorXd avgSNRs = Eigen::VectorXd(numOfTests);
	Eigen::VectorXd avgSARs = Eigen::VectorXd(numOfTests);

	for(int i=0; i<numOfTests; i++)
	{
		//Get samples from each sound
		Eigen::VectorXd ilSamples;
		Eigen::VectorXd irSamples;
		Eigen::VectorXd jlSamples;
		Eigen::VectorXd jrSamples;
		sampleFiles->data()[samples0[i]-1]->readSoundSample(&ilSamples, &irSamples);
		sampleFiles->data()[samples1[i]-1]->readSoundSample(&jlSamples, &jrSamples);
		
		long maxLength = max(ilSamples.size(), max(jlSamples.size(), noiselSamples.size()));
		Eigen::VectorXd zeroVector = Eigen::VectorXd::Zero(maxLength);
		noiselSamples.conservativeResizeLike(zeroVector);
		noiserSamples.conservativeResizeLike(zeroVector);
		ilSamples.conservativeResizeLike(zeroVector);
		irSamples.conservativeResizeLike(zeroVector);
		jlSamples.conservativeResizeLike(zeroVector);
		jrSamples.conservativeResizeLike(zeroVector);

		//Generate stereo positions
		double iStereo = 0.5 * (1 + stereoSep[i]);
		double jStereo = 0.5 * (1 - stereoSep[i]);
		//double iStereo = 1.;
		//double jStereo = 0.;

		//Pan sounds
		Eigen::VectorXd iSamples = (ilSamples + irSamples)/sqrt(2.);
		ilSamples = sqrt(1-iStereo) * iSamples;
		irSamples = sqrt(iStereo) * iSamples;
		Eigen::VectorXd jSamples = (jlSamples + jrSamples)/sqrt(2.);
		jlSamples = sqrt(1-jStereo) * jSamples;
		jrSamples = sqrt(jStereo) * jSamples;

		//Mix sounds
		Eigen::VectorXd mixlSamples = noiselSamples + ilSamples + jlSamples;
		Eigen::VectorXd mixrSamples = noiserSamples + irSamples + jrSamples;

		inSNRs(i) = 10*log10(((ilSamples + jlSamples).squaredNorm() + (irSamples + jrSamples).squaredNorm()) / (noiselSamples.squaredNorm() + noiserSamples.squaredNorm()));
		cout << "inSNR: " << inSNRs(i) << endl;

		//Perform separation
		srand(19873495);
		vector<vector<Eigen::VectorXd>>* separated = Separation::separate(&mixlSamples, &mixrSamples, 2, fOp, Separation::SOFT, false, false);
		long recLength = separated->data()[0][0].size();

		//Concatenate left and right channels
		Eigen::MatrixXd concatOriginals = Eigen::MatrixXd(2*recLength, 2);
		zeroVector = Eigen::VectorXd::Zero(recLength);
		ilSamples.conservativeResizeLike(zeroVector);
		irSamples.conservativeResizeLike(zeroVector);
		concatOriginals.col(0) << ilSamples, irSamples;
		jlSamples.conservativeResizeLike(zeroVector);
		jrSamples.conservativeResizeLike(zeroVector);
		concatOriginals.col(1) << jlSamples, jrSamples;

		for(int sep=0; sep<2; sep++)
		{
			for(int channel=0; channel<2; channel++)
			{
				for(int t=0; t<recLength; t++)
				{
					if(separated->data()[sep][channel](t) != separated->data()[sep][channel](t))
					{
						separated->data()[sep][channel](t) = 0.;
					}
				}
			}
		}

		//cout << separated->data()[0][1](recLength-1) << endl;

		Eigen::MatrixXd concatSeparated = Eigen::MatrixXd(2*recLength, 2);
		for(int sep=0; sep<2; sep++)
		{
			Eigen::VectorXd concat = Eigen::VectorXd(2*recLength);
			concat << separated->data()[sep][0], separated->data()[sep][1];
			concatSeparated.col(sep) = concat;
		}

		/*if (i == 13)
		{
			noiseFile.writeDerivedOutput(0, &(separated->data()[0][0]), &(separated->data()[0][1]));
			noiseFile.writeDerivedOutput(1, &(separated->data()[1][0]), &(separated->data()[1][1]));
		}*/

		delete separated;

		noiselSamples.conservativeResizeLike(zeroVector);
		noiserSamples.conservativeResizeLike(zeroVector);
		Eigen::VectorXd concatNoise = Eigen::VectorXd(2*recLength);
		concatNoise << noiselSamples, noiserSamples;

		Eigen::MatrixXd Rss = concatOriginals.transpose() * concatOriginals;
		Eigen::MatrixXd projections = concatSeparated.transpose() * concatOriginals;
		Eigen::VectorXd noiseProj = concatSeparated.transpose() * concatNoise;
		Eigen::MatrixXd expansions = projections * Rss.inverse().conjugate().transpose() * concatOriginals.transpose();
		Eigen::MatrixXd SDR = Eigen::MatrixXd(2, 2);
		Eigen::MatrixXd SIR = Eigen::MatrixXd(2, 2);
		Eigen::VectorXd SNR = Eigen::VectorXd(2);
		Eigen::MatrixXd SAR = Eigen::MatrixXd(2, 2);
		for(int out=0; out<2; out++)
		{
			for(int in=0; in<2; in++)
			{
				Eigen::VectorXd s_target = (projections(out, in)/Rss(in, in)) * concatOriginals.col(in);
				Eigen::VectorXd e_interf = expansions.row(out).transpose() - s_target;
				Eigen::VectorXd e_noise = (noiseProj(out) / concatNoise.squaredNorm()) * concatNoise;
				Eigen::VectorXd e_artef = concatSeparated.col(out) - expansions.row(out).transpose() - e_noise;

				SDR(out, in) = 10 * log10(s_target.squaredNorm() / (e_interf + e_noise + e_artef).squaredNorm());
				SIR(out, in) = 10 * log10(s_target.squaredNorm() / e_interf.squaredNorm());
				SNR(out) = 10 * log10((s_target + e_interf).squaredNorm() / e_noise.squaredNorm());
				SAR(out, in) = 10 * log10((s_target + e_interf + e_noise).squaredNorm() / e_artef.squaredNorm());
			}
		}

		avgSDRs(i) = max(SDR(0, 0) + SDR(1, 1), SDR(0, 1) + SDR(1, 0))/2;
		avgSIRs(i) = max(SIR(0, 0) + SIR(1, 1), SIR(0, 1) + SIR(1, 0))/2;
		avgSNRs(i) = (SNR(0) + SNR(1))/2;
		avgSARs(i) = max(SAR(0, 0) + SAR(1, 1), SAR(0, 1) + SAR(1, 0))/2;

		cout << "SNR: " << SNR << endl;
	}


	cout << endl;
	double meanInSNR = inSNRs.mean();
	double sigmaInSNR = sqrt((inSNRs.array() - meanInSNR).square().sum() / (numOfTests-1));
	cout << "InSNR: " << meanInSNR << "+-" << sigmaInSNR << endl;
	sort(inSNRs.data(), inSNRs.data() + inSNRs.size());
	cout << inSNRs(0) << ", " << inSNRs(numOfTests/4) << ", " << inSNRs(numOfTests/2) << ", " << inSNRs(3*numOfTests/4) << ", " << inSNRs(numOfTests-1) << endl;
	double meanSDR = avgSDRs.mean();
	double sigmaSDR = sqrt((avgSDRs.array() - meanSDR).square().sum() / (numOfTests-1));
	cout << "SDR: " << meanSDR << " +- " << sigmaSDR << endl;
	sort(avgSDRs.data(), avgSDRs.data() + avgSDRs.size());
	cout << avgSDRs(0) << ", " << avgSDRs(numOfTests/4) << ", " << avgSDRs(numOfTests/2) << ", " << avgSDRs(3*numOfTests/4) << ", " << avgSDRs(numOfTests-1) << endl;
	double meanSIR = avgSIRs.mean();
	double sigmaSIR = sqrt((avgSIRs.array() - meanSIR).square().sum() / (numOfTests-1));
	cout << "SIR: " << meanSIR << " +- " << sigmaSIR << endl;
	sort(avgSIRs.data(), avgSIRs.data() + avgSIRs.size());
	cout << avgSIRs(0) << ", " << avgSIRs(numOfTests/4) << ", " << avgSIRs(numOfTests/2) << ", " << avgSIRs(3*numOfTests/4) << ", " << avgSIRs(numOfTests-1) << endl;
	double meanSNR = avgSNRs.mean();
	double sigmaSNR = sqrt((avgSNRs.array() - meanSNR).square().sum() / (numOfTests-1));
	cout << "SNR: " << meanSNR << " +- " << sigmaSNR << endl;
	sort(avgSNRs.data(), avgSNRs.data() + avgSNRs.size());
	cout << avgSNRs(0) << ", " << avgSNRs(numOfTests/4) << ", " << avgSNRs(numOfTests/2) << ", " << avgSNRs(3*numOfTests/4) << ", " << avgSNRs(numOfTests-1) << endl;
	double meanSAR = avgSARs.mean();
	double sigmaSAR = sqrt((avgSARs.array() - meanSAR).square().sum() / (numOfTests-1));
	cout << "SAR: " << meanSAR << " +- " << sigmaSAR << endl;
	sort(avgSARs.data(), avgSARs.data() + avgSARs.size());
	cout << avgSARs(0) << ", " << avgSARs(numOfTests/4) << ", " << avgSARs(numOfTests/2) << ", " << avgSARs(3*numOfTests/4) << ", " << avgSARs(numOfTests-1) << endl;

	for(int i=0; i<sampleFiles->size(); i++)
	{
		delete sampleFiles->data()[i];
	}
	delete sampleFiles;
}
