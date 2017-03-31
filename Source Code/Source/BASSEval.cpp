#include "BASSEval.h"
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

	double sumSDR = 0;
	double sumSIR = 0;
	double sumSAR = 0;
	for(int i=0; i<numOfSounds; i++)
	{
		double maxSDR = numeric_limits<double>::lowest();
		int mapTo = -1;
		for(int j=0; j<numOfSounds; j++)
		{
			double sdr = SDR(i, j);
			if(maxSDR < sdr)
			{
				maxSDR = sdr;
				mapTo = j;
			}
		}

		sumSDR += SDR(i, mapTo);
		sumSIR += SIR(i, mapTo);
		sumSAR += SAR(i, mapTo);
	}

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