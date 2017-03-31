/**
* CLASS: WavFileManager
* AUTHOR: WVS22
*
* Each WavFileManager instance handles a single sound sample from C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Sound Samples/.
* The data samples can be read from these files and (after processing) derivative files can be written in C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Output Sounds SM/.
*/
#pragma once
#include "../JuceLibraryCode/JuceHeader.h"
#include <vector>
#include <Eigen/Core>
using namespace std;
class WavFileManager
{
	//WVS22 - 08/01/2017 - Get samples from file in chunks
	static const int samplesPerChunk = 1024;
	juce::AudioFormatManager formatManager;
	juce::String filename;
	juce::AudioFormatReader* reader;
	juce::StringPairArray metadata;
public:
	WavFileManager(juce::String iFilename);
	~WavFileManager(void);

	/**
	* WavFileManager::readSoundSample(vector<double>&, vector<double>&)
	*
	* Reads the contents of the wav file and pushes them to the end of the vectors passed in.
	* Assumes that the file has 2 channels, left and right (mono channel is duplicated).
	* NOTE: least significant difference between samples will be 65536, corresponding to 
	* 3.051757812500000e-05 in MATLAB. Consider this scaling when comparing values to MATLAB.
	*/
	void readSoundSample(vector<double>& lSampleVector, vector<double>& rSampleVector);
	void readSoundSample(Eigen::VectorXd* lSampleVector, Eigen::VectorXd* rSampleVector);

	/**
	* WavFileManager::writeDerivedOutput(int, vector<double>&, vector<double>&)
	*
	* Creates a new wav file in C:/Users/Will/OneDrive/Uni/Individual Project/SoundSeparation/Output Sounds SM/
	* entitled sampleName_outputNumber.wav containing the samples provided.
	* Assumes that the sample vectors have the same length (uses left vector to determine
	* output length). Uses the same sample rate, bits per sample and metadata as the input
	* file.
	*/
	void writeDerivedOutput(int outputNumber, vector<double>& lSampleVector, vector<double>& rSampleVector);
	void writeDerivedOutput(int outputNumber, Eigen::VectorXd* lSampleVector, Eigen::VectorXd* rSampleVector);

};
