#include "Transform.h"
#include <cmath>

#define PI 3.14159265358979323846

const double Transform::nmfStepThreshold = 1e-7;

double* Transform::hamming(long size)
{
	//Create array for resule
	double* retVal = new double[size];

	//Calculate samples of periodic Hamming window
	for(long i=0; i<size; i++)
	{
		retVal[i] = 0.54 - 0.46*cos(2*PI*i/size);
	}

	return retVal;
}

//Helper function for Transform::FFT
//Recursively compute FFT in-place using butterfly network
//Assumes vector is of length 2^n
void FFTstep(complex<double>* vector, int n)
{
	//Handle base case
	if (n == 0) return;

	//Get parameters
	const int N = 1<<n;
	const complex<double> w = complex<double>(cos(2*PI/N), -sin(2*PI/N));

	//Perform butterfly mixing
	for (int i=0; i<N/2; i++)
	{
		complex<double> a = vector[i];
		complex<double> b = vector[i+N/2];
		vector[i] = a+b;
		vector[i+N/2] = a-b;
	}

	//Multiply by exponential
	complex<double> W = w;
	for (int k=N/2+1; k<N; k++)
	{
		vector[k] = vector[k] * W;
		W = W * w;
	}
	
	//Recurse on each half
	FFTstep(vector, n-1);
	FFTstep(vector+(N/2), n-1);
}

//Helper function for Transform::FFT
//Reverses the lower bits-many bits from value
//e.g. revbits(5, 5) => 20
long revbits(long value, int bits)
{
	long result = 0;
	for (int i=0; i<bits; i++)
	{
		result *= 2;
		result += value%2;
		value /= 2;
	}
	return result;
}

vector<complex<double>>* Transform::FFT(complex<double>* samples, long length)
{
	//Get size parameters
	int n = (int) (log(length-1.)/log(2.)) + 1;
	long paddedLength = 1<<n;

	//Pad samples
	vector<complex<double>>* retVal = new vector<complex<double>>(paddedLength);
	for(long i=0; i<length; i++)
	{
		(*retVal)[i] = samples[i];
	}

	//Perform in-place butterfly FFT
	FFTstep(retVal->data(), n);

	//Re-arrange butterfly outputs
	for (long i=0; i<paddedLength; i++)
	{
		long revi = revbits(i, n);
		if (i<revi)
		{
			complex<double> temp = (*retVal)[i];
			(*retVal)[i] = (*retVal)[revi];
			(*retVal)[revi] = temp;
		}
	}

	return retVal;
}

vector<complex<double>>* Transform::halfFFT(double* samples, long nfft)
{
	//Since input is real, can pack into a half-size FFT
	const long packedLength = (nfft+1)/2;
	complex<double>* packedSamples = new complex<double>[packedLength];
	for(long i=0; i<nfft/2; i++)
	{
		packedSamples[i] = complex<double>(samples[2*i], samples[2*i + 1]);
	}
	if(packedLength > nfft/2)
	{
		packedSamples[packedLength-1] = complex<double>(samples[nfft-1]);
	}

	//Pad and compute FFT
	vector<complex<double>>* fft = FFT(packedSamples, packedLength);
	long fftLength = fft->size();

	//Unpack values up to Nyquist limit
	vector<complex<double>>* retVal = new vector<complex<double>>(fftLength+1);
	(*retVal)[0] = ((*fft)[0] + conj((*fft)[0]) - complex<double>(0.,1.)*((*fft)[0]-conj((*fft)[0])))*0.5;
	const complex<double> w = complex<double>(cos(PI/fftLength), -sin(PI/fftLength));
	complex<double> W = w;
	for(long i=1; i<fftLength; i++)
	{
		complex<double> Fe = (*fft)[i] + conj((*fft)[fftLength-i]);
		complex<double> Fo = (*fft)[i] - conj((*fft)[fftLength-i]);
		(*retVal)[i] = (Fe - complex<double>(0.,1.)*W*Fo) * 0.5;
		W = W*w;
	}
	(*retVal)[fftLength] = ((*fft)[0] + conj((*fft)[0]) + complex<double>(0.,1.)*((*fft)[0] - conj((*fft)[0])))*0.5;

	//Clean up
	delete packedSamples;
	delete fft;

	return retVal;
}

double* Transform::ifft(vector<complex<double>>* fft)
{
	//Inverse DFT can be done by applying conjugate, then forwards DFT, conjugate and scale
	//Conjugate input
	const long length = fft->size();
	complex<double>* fftConj = new complex<double>[length];
	for(long i=0; i<length; i++)
	{
		fftConj[i] = conj((*fft)[i]);
	}

	//Compute FFT forwards
	vector<complex<double>>* complexSamples = FFT(fftConj, length);

	//Since only considering real output, conjugate does nothing, so just take the real part
	double* retVal = new double[length];
	for(long i=0; i<length; i++)
	{
		retVal[i] = (*complexSamples)[i].real()/length;
	}

	//Clean up
	delete fftConj;
	delete complexSamples;

	return retVal;
}


Matrix<complex<double>>* Transform::stft(vector<double>* samples, long windowSize, long hopSize, long nfft)
{
	//Get spectrum size
	const long sampleCount = samples->size();
	const long columns = 1 + (sampleCount - windowSize)/hopSize;

	//Get access to samples and window function
	double* sampleP = samples->data();
	double* window = hamming(windowSize);

	//Calculate spectrum
	Matrix<complex<double>>* spectrum = new Matrix<complex<double>>(0, 0);
	double* windowedFrame = new double[nfft];
	for(long col=0; col<columns; col++)
	{
		//In each time frame, find the start and apply window function
		const long startOfFrame = col*hopSize;
		double* frame = sampleP + startOfFrame;
		for(long i=0; i<windowSize; i++)
		{
			windowedFrame[i] = frame[i] * window[i];
		}

		//Apply FFT to get spectrum
		//Only need half FFT as data is real
		vector<complex<double>>* columnData = halfFFT(windowedFrame, nfft);
		spectrum->data.push_back(*columnData);
		delete columnData;
	}

	//Clean up
	delete window;
	delete windowedFrame;

	return spectrum;
}

vector<double>* Transform::istft(Matrix<complex<double>>* stft, long hopSize)
{
	//Get parameters
	const long nfft = stft->data[0].size();
	const long columns = stft->data.size();
	const long sampleCount = nfft + (columns-1)*hopSize;

	//Create vector for output samples
	vector<double>* retVal = new vector<double>(sampleCount);

	//Get window function
	double* window = hamming(nfft);

	//Perform ISTFT
	for(long col=0; col<columns; col++)
	{
		//For each time step, find the position to place the samples at
		const long startOfFrame = col*hopSize;

		//Calculate the samples for the time frame
		//Can use the half FFT directly since we are only taking the real part of the result
		double* frame = ifft(&(stft->data[col]));

		//Apply window function to samples to interpolate between frames and add to result
		for(long i=0; i<nfft; i++)
		{
			(*retVal)[startOfFrame+i] += frame[i] * window[i];
		}

		delete frame;
	}

	//Remove scaling caused by window function
	double sum = 0.;
	for(long i=0; i<nfft; i++)
	{
		sum += window[i]*window[i];
	}
	const double scale = hopSize/sum;
	for(long i=0; i<sampleCount; i++)
	{
		(*retVal)[i] *= scale;
	}

	//Clean up
	delete window;

	return retVal;
}

vector<Matrix<double>>* Transform::nmf(Matrix<double>* x, long n)
{
	//Get matrix parameters
	long tCount = x->data.size();
	long fCount = x->data[0].size();

	//Select random rows and columns from x for initial mix and source estimates
	Matrix<double> mix = Matrix<double>(tCount, n);
	Matrix<double> source = Matrix<double>(0, fCount);
	for(long i=0; i<n; i++)
	{
		long randIndex = rand() % tCount;
		source.data.push_back(vector<double>(x->data[randIndex]));
		randIndex = rand() % fCount;
		for(long t=0; t<tCount; t++)
		{
			mix.data[t][i] = x->data[t][randIndex];
		}
	}

	//Identify scale of x
	double xScale = 0.;
	for(long t=0; t<tCount; t++)
	{
		for(long f=0; f<fCount; f++)
		{
			xScale += x->data[t][f];
		}
	}

	//Iteratively improve estimates
	double lastDistance = numeric_limits<double>::max();
	while(true)
	{
		//Compute matrix multiplications for multiplicative updates
		Matrix<double>* mixSource = Matrix<double>::multiply(&mix, &source);
		source.transpose();
		Matrix<double>* mixModNumerator = Matrix<double>::multiply(x, &source);
		Matrix<double>* mixModDenominator = Matrix<double>::multiply(mixSource, &source);
		source.transpose();
		mix.transpose();
		Matrix<double>* sourceModNumerator = Matrix<double>::multiply(&mix, x);
		Matrix<double>* sourceModDenominator = Matrix<double>::multiply(&mix, mixSource);
		mix.transpose();

		//Perform multiplicative updates
		for(long t=0; t<tCount; t++)
		{
			for(long i=0; i<n; i++)
			{
				if (mixModDenominator->data[t][i] == 0)
				{
					mixModDenominator->data[t][i] = numeric_limits<double>::lowest();
				}
				mix.data[t][i] *= mixModNumerator->data[t][i] / mixModDenominator->data[t][i];
			}
		}
		for(long i=0; i<n; i++)
		{
			for(long f=0; f<fCount; f++)
			{
				if (sourceModDenominator->data[i][f] == 0)
				{
					sourceModDenominator->data[i][f] = 1e-307;
				}
				source.data[i][f] *= sourceModNumerator->data[i][f] / sourceModDenominator->data[i][f];
			}
		}

		//Clean up intermediate matrices
		delete mixSource;
		delete mixModNumerator;
		delete mixModDenominator;
		delete sourceModNumerator;
		delete sourceModDenominator;

		//Identify scales of current estimate of x
		Matrix<double>* currentEstimate = Matrix<double>::multiply(&mix, &source);
		double eScale = 0.;
		for(long t=0; t<tCount; t++)
		{
			for(long f=0; f<fCount; f++)
			{
				eScale += currentEstimate->data[t][f];
			}
		}

		//Calculate Euclidean squared distance between x and estimate, considering scales
		double distance = 0.;
		for(long t=0; t<tCount; t++)
		{
			for(long f=0; f<fCount; f++)
			{
				distance += pow((x->data[t][f] / xScale) - (currentEstimate->data[t][f] / eScale), 2);
			}
			for(long i=0; i<n; i++)
			{
				mix.data[t][i] *= xScale / eScale;
			}
		}
		delete currentEstimate;

		//Stop iterating if making little progress
		if(lastDistance - distance < nmfStepThreshold)
		{
			break;
		}
		lastDistance = distance;
	}

	//Compile return value
	vector<Matrix<double>>* retVal = new vector<Matrix<double>>();
	retVal->push_back(mix);
	retVal->push_back(source);

	return retVal;
}
