#include "WavFileManager.h"
#include <iostream>

WavFileManager::WavFileManager(juce::String iFilename)
{
	formatManager.registerBasicFormats();
	//sampleName = iSampleName;
	//juce::String filename = "C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/" + sampleName + ".wav";
	filename = iFilename;
	juce::File inFile(filename);
	reader = formatManager.createReaderFor(inFile);
	metadata = reader->metadataValues;
}


WavFileManager::~WavFileManager(void)
{
	delete reader;
}

void WavFileManager::readSoundSample(vector<double>& lSampleVector, vector<double>& rSampleVector)
{
	//Get parameters
	const long sampleCount = reader->lengthInSamples;
	const int numberOfChunks = sampleCount/samplesPerChunk;
	
	//Create space for chunk buffer
	int **sampleChunk = new int*[2];
	sampleChunk[0] = new int[samplesPerChunk];
	sampleChunk[1] = new int[samplesPerChunk];

	//Obtain each chunk
	for(int i=0; i<numberOfChunks; i++)
	{
		int startSample = i*samplesPerChunk;
		reader->read(sampleChunk,2,startSample,samplesPerChunk,true);
		for(int j=0; j<samplesPerChunk; j++)
		{
			lSampleVector.push_back((double)sampleChunk[0][j]);
			rSampleVector.push_back((double)sampleChunk[1][j]);
		}
	}

	//Process any remaining samples
	const int lastChunkStart = numberOfChunks*samplesPerChunk;
	const int remainder = sampleCount - lastChunkStart;
	reader->read(sampleChunk,2,lastChunkStart,remainder,true);
	for(int j=0; j<remainder; j++)
	{
		lSampleVector.push_back((double)sampleChunk[0][j]);
		rSampleVector.push_back((double)sampleChunk[1][j]);
	}

	//Clean up
	delete sampleChunk[0];
	delete sampleChunk[1];
	delete sampleChunk;
}

void WavFileManager::readSoundSample(Eigen::VectorXd* lSampleVector, Eigen::VectorXd* rSampleVector)
{
	vector<double> lSamples = vector<double>();
	vector<double> rSamples = vector<double>();

	//Get parameters
	const long sampleCount = reader->lengthInSamples;
	const int numberOfChunks = sampleCount/samplesPerChunk;
	
	//Create space for chunk buffer
	int **sampleChunk = new int*[2];
	sampleChunk[0] = new int[samplesPerChunk];
	sampleChunk[1] = new int[samplesPerChunk];

	//Obtain each chunk
	for(int i=0; i<numberOfChunks; i++)
	{
		int startSample = i*samplesPerChunk;
		reader->read(sampleChunk,2,startSample,samplesPerChunk,true);
		for(int j=0; j<samplesPerChunk; j++)
		{
			lSamples.push_back((double)sampleChunk[0][j]);
			rSamples.push_back((double)sampleChunk[1][j]);
		}
	}

	//Process any remaining samples
	const int lastChunkStart = numberOfChunks*samplesPerChunk;
	const int remainder = sampleCount - lastChunkStart;
	reader->read(sampleChunk,2,lastChunkStart,remainder,true);
	for(int j=0; j<remainder; j++)
	{
		lSamples.push_back((double)sampleChunk[0][j]);
		rSamples.push_back((double)sampleChunk[1][j]);
	}

	(*lSampleVector) = Eigen::VectorXd::Map(lSamples.data(), lSamples.size());
	(*rSampleVector) = Eigen::VectorXd::Map(rSamples.data(), rSamples.size());

	//Clean up
	delete sampleChunk[0];
	delete sampleChunk[1];
	delete sampleChunk;
}

//WVS22 - 08/01/2017 - Saving samples to output
void WavFileManager::writeDerivedOutput(int outputNumber, vector<double>& lSampleVector, vector<double>& rSampleVector)
{
	//Get writer for the file
	juce::WavAudioFormat wavFormat;
	//juce::String filename = "C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Output Sounds SM/" + sampleName + "_" + juce::String(outputNumber) + ".wav";
	juce::String outFilename = filename.dropLastCharacters(4) + "_" + juce::String(outputNumber) + ".wav";
	juce::File outFile(outFilename);
	outFile.deleteFile();
	juce::FileOutputStream* outputStream = outFile.createOutputStream();
	const double sampleRate = reader->sampleRate;
	const int bitsPerSample = reader->bitsPerSample;
	juce::AudioFormatWriter* writer = wavFormat.createWriterFor(outputStream,sampleRate,2,bitsPerSample,metadata,0);

	//Format samples
	const int sampleCount = lSampleVector.size();
	int** intSamples = new int*[2];
	intSamples[0] = new int[sampleCount];
	intSamples[1] = new int[sampleCount];
	for(int i=0; i<sampleCount; i++)
	{
		intSamples[0][i] = (int)lSampleVector[i];
		intSamples[1][i] = (int)rSampleVector[i];
	}

	//Write samples to file
	const int** constIntSamples = new const int*[2];
	constIntSamples[0] = intSamples[0];
	constIntSamples[1] = intSamples[1];
	writer->write(constIntSamples,sampleCount);
	writer->flush();

	//Clean up
	delete intSamples[0];
	delete intSamples[1];
	delete intSamples;
	delete constIntSamples;
	delete writer;
}

void WavFileManager::writeDerivedOutput(int outputNumber, Eigen::VectorXd* lSampleVector, Eigen::VectorXd* rSampleVector)
{
	//Get writer for the file
	juce::WavAudioFormat wavFormat;
	//juce::String filename = "C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Output Sounds SM/" + sampleName + "_" + juce::String(outputNumber) + ".wav";
	juce::String outFilename = filename.dropLastCharacters(4) + "_" + juce::String(outputNumber) + ".wav";
	juce::File outFile(outFilename);
	outFile.deleteFile();
	juce::FileOutputStream* outputStream = outFile.createOutputStream();
	const double sampleRate = reader->sampleRate;
	const int bitsPerSample = reader->bitsPerSample;
	juce::AudioFormatWriter* writer = wavFormat.createWriterFor(outputStream,sampleRate,2,bitsPerSample,metadata,0);

	//Format samples
	const int sampleCount = lSampleVector->size();
	int** intSamples = new int*[2];
	intSamples[0] = new int[sampleCount];
	intSamples[1] = new int[sampleCount];
	for(int i=0; i<sampleCount; i++)
	{
		intSamples[0][i] = (int)(*lSampleVector)(i);
		intSamples[1][i] = (int)(*rSampleVector)(i);
	}

	//Write samples to file
	const int** constIntSamples = new const int*[2];
	constIntSamples[0] = intSamples[0];
	constIntSamples[1] = intSamples[1];
	writer->write(constIntSamples,sampleCount);
	writer->flush();

	//Clean up
	delete intSamples[0];
	delete intSamples[1];
	delete intSamples;
	delete constIntSamples;
	delete writer;
}
