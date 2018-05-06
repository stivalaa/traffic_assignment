#include "microblas_double.h"

#include <math.h>

/* -------------------  GENERALIZED SPARSE DDOT ------------------------ */
/*                                                                       */
/*                          r = dot(x,y);                                */
/*                                                                       */
/* ----------------------------------------------------------------------*/

double ublas_spdotab_double(int nz,
        const double *val, const int *ind, const double *y, int incy)
{ 
    double t= 0.0;
    int i=0; 
    for (; i<nz; i++) 
        t += val[i] * y[ind[i]*incy]; 
    return t;
}

/* -----------------  DIAGONAL SPARSE MULTPLICATION -------------------- */
/*                                                                       */
/*                      y = alpha * diag(A) * x  + y                     */
/*                                                                       */
/* ----------------------------------------------------------------------*/
/*  diags is vector, int(M), that point to where the diagonals are in    */
/*  val[] array.                                                         */
/* ----------------------------------------------------------------------*/


void ublas_diagmult_double(int M, double alpha, const double *diags,  
            const double *x, int incx, double *y, int incy)
{

    int i;

    for (i=0; i<M; i++)
      y[i*incy] =  fmin( y[i*incy],  alpha + diags[i] + x[i*incx] ); /* (min,+) algebra */
/*        y[i*incy] +=  alpha * diags[i] * x[i*incx]; */
}



/* ----------------------  GENERALIZED DAPXBY -------------------------- */
/*                                                                       */
/*                           y = alpha*x + y                             */
/*                                                                       */
/* ----------------------------------------------------------------------*/

void ublas_spaxpy_double(double alpha, int nz, 
    const double *val, const int *ind, double *y, int incy) 
{ 
    int i=0; 
    for ( ; i<nz; i++) 
            y[ind[i] * incy] += alpha * val[i]; 
}

    
