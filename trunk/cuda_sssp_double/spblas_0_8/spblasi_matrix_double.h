
#ifndef SPBLAS_MATRIX_double_H
#define SPBLAS_MATRIX_double_H

#include "blas_enum.h"
#include "spblasi_matrix.h"
#include "csr_double.h"

/**
	An SPBLASI_Matrix_double is just a generic SPBLASI_Matrix strucutre whose
	components have already been identified as double-precision reals.  This
	allows use to more precisely describe the argument when being passed into
	functions.
*/
typedef SPBLASI_Matrix SPBLASI_Matrix_double;

/*
* SPBLASI_Matrix_double defines the combined structure to support
* the various formats of sparse matrices used in this implementation
* of the Sparse BLAS.
*
*/

#define SPBLASI_Matrix_double_diag(A)  ((double *) SPBLASI_Matrix_diag(A))
#define SPBLASI_Matrix_double_CSR(A)  ((CSR_double *) SPBLASI_Matrix_CSR(A))

											
SPBLASI_Matrix *SPBLASI_Matrix_double_new(int M, int N, 
					SPBLASI_Matrix_prop_type prop); 

SPBLASI_Matrix *SPBLASI_Matrix_double_new_block(int Mb, int Nb, int k, int l, 
						SPBLASI_Matrix_prop_type prop);

SPBLASI_Matrix *SPBLASI_Matrix_double_new_variable_block(int Mb, int Nb, 
		const int *K, const int *L, SPBLASI_Matrix_prop_type prop);

int SPBLASI_Matrix_double_prepare_for_insert(SPBLASI_Matrix *A);

int SPBLASI_Matrix_double_insert_entry( SPBLASI_Matrix *A, double val, 
					int i, int j);

int SPBLASI_Matrix_double_insert_entries(SPBLASI_Matrix *A, 
		int nz, const double *Val, 
		const int *I, const int *J);

int SPBLASI_Matrix_double_insert_col( SPBLASI_Matrix *A, 
		int col, int nz, const double *Val, 
						const int *I);

int SPBLASI_Matrix_double_insert_row( SPBLASI_Matrix *A, 
		int row, int nz, const double *Val, 
						const int *J);

int SPBLASI_Matrix_double_insert_block(SPBLASI_Matrix *A, 
		const double *Val, int row_stride, 
					int col_stride, int bi, int bj);

int SPBLASI_Matrix_double_insert_clique( SPBLASI_Matrix *A, 
		int k, int l, const double *Val, 
				int row_stride, int col_stride, const int *I, const int *J);

int SPBLASI_Matrix_double_end_construction(SPBLASI_Matrix *A);

int SPBLASI_Matrix_double_delete(SPBLASI_Matrix *A);

int SPBLASI_Matrix_double_usmm(enum blas_order_type order, 
			enum blas_trans_type trans, 
			int k, double alpha, const SPBLASI_Matrix *A, 
			const double *B, int lbB, 
			const double *C, int ldC);

int SPBLASI_Matrix_double_usmv(
			enum blas_trans_type trans, 
			double alpha, const SPBLASI_Matrix *A, 
			const double *x, int incx, 
			double *y, int incy);

int SPBLASI_Matrix_double_ussm(enum blas_order_type order, 
			enum blas_trans_type trans, int k, double alpha, 
			const SPBLASI_Matrix *A, 
			double *B, int ldB);

int SPBLASI_Matrix_double_ussv( enum blas_trans_type trans, 
  			double alpha, 
			const SPBLASI_Matrix *A, 
			double *x, int incx);


/* Diagnostic tools */
void SPBLASI_Matrix_double_print(const SPBLASI_Matrix *A, const char *format);

#endif
