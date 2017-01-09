#include "../JuceLibraryCode/JuceHeader.h"

//WVS22 - 03/01/2017 - Include Headers
#include <iostream>
#include <string>

//WVS22 - 08/01/2017 - Getting samples from input
#include <vector>

using namespace std;

//WVS22 - 08/01/2017 - Get samples from file in chunks
static const int chunkSize = 1024;

//WVS22 - 08/01/2017 - Getting samples from input
void readSamplesFromReader(juce::AudioFormatReader* reader, vector<double>& lSampleVector, vector<double>& rSampleVector)
{
	//Get parameters
	const int sampleCount = reader->lengthInSamples;
	const int numberOfChunks = sampleCount/chunkSize;

	//Create space for chunk buffer
	int **sampleChunk = new int*[2];
	sampleChunk[0] = new int[chunkSize];
	sampleChunk[1] = new int[chunkSize];

	//Obtain each chunk
	for(int i=0; i<numberOfChunks; i++)
	{
		int startSample = i*chunkSize;
		reader->read(sampleChunk,2,startSample,chunkSize,true);
		for(int j=0; j<chunkSize; j++)
		{
			lSampleVector.push_back((double)sampleChunk[0][j]);
			rSampleVector.push_back((double)sampleChunk[1][j]);
		}
	}

	//Process any remaining samples
	const int lastChunkStart = numberOfChunks*chunkSize;
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

//WVS22 - 08/01/2017 - Saving samples to output
void writeSamplesToWriter(juce::AudioFormatWriter* writer, vector<double>& lSampleVector, vector<double>& rSampleVector)
{
	const int sampleCount = lSampleVector.size();
	int** intSamples = new int*[2];
	intSamples[0] = new int[sampleCount];
	intSamples[1] = new int[sampleCount];
	for(int i=0; i<sampleCount; i++)
	{
		intSamples[0][i] = (int)lSampleVector[i];
		intSamples[1][i] = (int)rSampleVector[i];
	}
	const int** constIntSamples = new const int*[2];
	constIntSamples[0] = intSamples[0];
	constIntSamples[1] = intSamples[1];
	writer->write(constIntSamples,sampleCount);
	writer->flush();

	delete intSamples[0];
	delete intSamples[1];
	delete intSamples;
}

int main (int argc, char* argv[])
{
	//WVS22 - 03/01/2017 - Test file input
	juce::AudioFormatManager formatManager;
	formatManager.registerBasicFormats();
	juce::File inFile("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/Guitar Trumpet/GD3vlfn_TGs41fn.wav");
	juce::AudioFormatReader* inReader = formatManager.createReaderFor(inFile);
	juce::StringPairArray metadata = inReader->metadataValues;

	cout << inReader->numChannels << " channels" << endl;
	cout << inReader->lengthInSamples << " samples" << endl;
	cout << inReader->sampleRate << " Hz" << endl;
	cout << inReader->bitsPerSample << " bits per sample" << endl;

	//false
	cout << inReader->usesFloatingPointData << endl;

	//WVS22 - 08/01/2017 - Getting samples from input
	const int sampleCount = inReader->lengthInSamples;
	vector<double> lSampleVector;
	vector<double> rSampleVector;
	readSamplesFromReader(inReader, lSampleVector, rSampleVector);

	cout << "Samples read" << endl;

	//WVS22 - 08/01/2017 - Saving samples to output
	juce::WavAudioFormat wavFormat;
	juce::File outFile("C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Output Sounds SM/Guitar Trumpet/GD3vlfn_TGs41fn.wav");
	outFile.deleteFile();
	juce::FileOutputStream* outputStream = outFile.createOutputStream();
	const int sampleRate = inReader->sampleRate;
	const int bitsPerSample = inReader->bitsPerSample;
	juce::AudioFormatWriter* outWriter = wavFormat.createWriterFor(outputStream,sampleRate,2,bitsPerSample,metadata,0);
	writeSamplesToWriter(outWriter, lSampleVector, rSampleVector);

	cout << "Samples written" << endl;

	string wait;
	cin >> wait;

	delete inReader;
	delete outWriter;

    return 0;
}


