#include <stdlib.h>
#include "blas_sparse.h"
#include "spblasi_error.h"  

/*
* Level 1 routines
*
*/

void BLAS_dusdot( enum blas_conj_type conj, int nz,
    const double *x,  const int *index,  const double *y, int incy,
    double *r, enum blas_base_type index_base)
{

	int i=0;
	double t= 0.0;

	if (index_base == blas_one_base)
		y -= incy;

	for (i=0; i<nz; i++)
		t += x[i] * y[index[i]*incy];

	*r = t;
}


		

void BLAS_dusaxpy(int nz, double alpha, const double *x,
    const int *index, double *y, int incy,
    enum blas_base_type index_base)
{
	int i=0;

	if (index_base == blas_one_base)
		y -= incy;

	for (i=0; i<nz; i++)
		y[index[i]*incy] += alpha * x[i];
}

void BLAS_dusga( int nz, const double *y, int incy, double *x, const int *indx,
              enum blas_base_type index_base )
{
	int i=0;

	if (index_base == blas_one_base)
		y -= incy;

	for (i=0; i<nz; i++)
		x[i] = y[indx[i]*incy];

}

void BLAS_dusgz( int nz, double *y, int incy, double *x, const int *indx,
              enum blas_base_type index_base )
{
	int i=0;

	if (index_base == blas_one_base)
		y -= incy;

	for (i=0; i<nz; i++)
	{
		x[i] = y[indx[i]*incy];
		y[indx[i]*incy] = 0.0;
	}

}


void BLAS_dussc(int nz, const double *x, double *y, int incy, const int *index,
    enum blas_base_type index_base)
{
	int i=0;

	if (index_base == blas_one_base)
		y -= incy;

	for (i=0; i<nz; i++)
		y[index[i]*incy] = x[i];

}



