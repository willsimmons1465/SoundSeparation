/**
* CLASS: Matrix
* AUTHOR: WVS22
*
* Instances represent two-dimensional matrices implemented using nested vectors.
*/
#pragma once
#include <vector>
using namespace std;
template <class T>
class Matrix
{
	bool transposed;
public:
	vector<vector<T>> data;

	/**
	* Matrix<T>::Matrix(long, long)
	*
	* Constructor for a new matrix. Initialises data to a dim0 x dim1 matrix.
	* Dimensions of data can still be changed by accessing it directly.
	*/
	Matrix<T>(long dim0, long dim1);

	/**
	* Matrix<T>::transpose(void)
	*
	* Transposes the matrix with regards to its behaviour in multiplication.
	* This does not affect the underlying data structure, so accessing the elements
	* will be the same.
	*/
	inline void transpose(void) {transposed = !transposed;}

	/**
	* Matrix<T>::multiply(Matrix<T>*, Matrix<T>*)
	*
	* Creates a new matrix corresponding to the result of multiplying the matrices
	* lhs and rhs. The transposed field is considered here.
	* Assumes that the supplied matrices match with the inner dimension of the
	* multiplication.
	*/
	static Matrix<T>* multiply(Matrix<T>* lhs, Matrix<T>* rhs);
};


template <class T>
Matrix<T>::Matrix(long dim0, long dim1) : transposed(false)
{
	data = vector<vector<T>>();
	for(long i=0; i<dim0; i++)
	{
		data.push_back(vector<T>(dim1));
	}
}

template <class T>
Matrix<T>* Matrix<T>::multiply(Matrix<T>* lhs, Matrix<T>* rhs)
{
	//Get dimension sizes
	long dim0 = (lhs->transposed)?lhs->data[0].size():lhs->data.size();
	long dim1 = (rhs->transposed)?rhs->data.size():rhs->data[0].size();
	long dimMult = (lhs->transposed)?lhs->data.size():lhs->data[0].size();

	//Naively multiply
	Matrix<T>* retVal = new Matrix<T>(dim0, dim1);
	for(long i0=0; i0<dim0; i0++)
	{
		for(long i1=0; i1<dim1; i1++)
		{
			T val = 0.;
			for(long iMult=0; iMult<dimMult; iMult++)
			{
				double lhsTerm = (lhs->transposed)?lhs->data[iMult][i0]:lhs->data[i0][iMult];
				double rhsTerm = (rhs->transposed)?rhs->data[i1][iMult]:rhs->data[iMult][i1];
				val += lhsTerm*rhsTerm;
			}
			retVal->data[i0][i1] = val;
		}
	}

	return retVal;
}
