/**
* CLASS: Transform
* AUTHOR: WVS22
*
* Defines a collection of static methods for mathematical functions.
* Typically operates on vectors of numbers when the length of the input is needed and simple arrays otherwise.
*/
#pragma once
#include <complex>
#include <vector>
#include "Matrix.h"
using namespace std;
class Transform
{
	static const double nmfStepThreshold;
public:

	/**
	* Transform::stft(vector<double>*, long, long, long)
	* 
	* Calculates the Short-Time-Fourier-Transform of the given sample vector according to
	* the given parameters. Windows using a Hamming window function.
	* In the output, the first dimension gives time and the second gives frequency bin,
	* e.g. (*s)[3][6] corresponds to the quantity in the 7th frequency bin for the 4th window.
	*/
	static Matrix<complex<double>>* stft(vector<double>* samples, long windowSize, long hopSize, long nfft);

	/**
	* Transform::istft(Matrix<complex<double>>*, long)
	*
	* Calculates the Inverse-Short-Time-Fourier-Transform (inverse of Transform::stft).
	* Assumes stft provided corresponds to a real input (so we are only given the half FFTs 
	* and can ignore the imaginary part of the result).
	*/
	static vector<double>* istft(Matrix<complex<double>>* stft, long hopSize);

	/**
	* Transform::hamming(long)
	*
	* Gives an array of length size containing the samples of a periodic Hamming window
	* of period size.
	*/
	static double* hamming(long size);

	/**
	* Transform::FFT(complex<double>*, long)
	*
	* Performs the Discrete Fourier Transform using FFT.
	* Assumes complex vector samples has size length.
	* Will automatically pad to next power of 2 >= length for efficient computation
	* so return value will likely be longer than length.
	*/
	static vector<complex<double>>* FFT(complex<double>* samples, long length);

	/**
	* Transform::halfFFT(double*, long)
	*
	* Performs the Discrete Fourier Transform using FFT and discards the upper half.
	* DFTs of real inputs have Hermitian symmetry so can be described entirely by either
	* half, including both the DC and Nyquist limit frequencies.
	* Assumes double vector samples has size nfft.
	* Uses Transform::FFT so will automatically pad up to the next power of 2 in size,
	* e.g. if nfft==11, return value will have size 9 (11->16, 16/2+1=9).
	*/
	static vector<complex<double>>* halfFFT(double* samples, long nfft);

	/**
	* Transform::ifft(vector<complex<double>>*)
	*
	* Performs the Inverse DFT using FFT.
	* Uses Transform::FFT so will automatically pad up to the next power of 2 in size.
	* Outputs only the real part of the time domain.
	* Designed to be inverse of Transform::halfFFT.
	*/
	static double* ifft(vector<complex<double>>* fft);

	/**
	* Transform::nmf(Matrix<double>*, long)
	*
	* Uses multiplicative-step Gradient descent to perform Non-negative Matrix Factorisation.
	* nmfStepThreshold describes the minimum improvement to Euclidean distance needed to
	* continue the descent; increase for more accuracy and decrease for speed.
	* Returns a pair of matrices, the first is the mixing matrix and the second is the source.
	*/
	static vector<Matrix<double>>* nmf(Matrix<double>* x, long factors);

};
