
#ifndef MICROBLAS_double_H
#define MICROBLAS_double_H


double ublas_spdotab_double(int nz, const double *val, const int *ind, 
		const double *y, int incy);

void ublas_spaxpy_double(double alpha, int nz,
        const double *val, const int *ind, 
        double *y, int incy) ;

void ublas_diagmult_double(int M, double alpha, const double *val, 
            const double *x, int incx, double *y, int incy);


#endif
