#include "Matrix.h"

/*template <class T>
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
	long dim0 = (lhs->transposed)?lhs->data[0].size():lhs->data.size();
	long dim1 = (rhs->transposed)?rhs->data.size():rhs->data[0].size();
	long dimMult = (lhs->transposed)?lhs->data.size():lhs->data[0].size();

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
}*/
