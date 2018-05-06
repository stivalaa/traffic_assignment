
#ifndef CSR_double_H
#define CSR_double_H 1

#include "spvec_double.h"

/**
  Compressed Row Sparse Matrix: based on list of sparse vectors 
  Modified sparse row stores diagonal entries in separate array.
*/

typedef struct
{
    int M_, N_, nz_;
	SPVEC_double **row_;		/* M sparse vectors: one for each row */


} CSR_double;

#define CSR_double_row(A,i)   ((A)->row_[i])
#define CSR_double_val(A, i, k)    SPVEC_double_val(CSR_double_row(A,i), k)
#define CSR_double_index(A, i, k)  SPVEC_double_index(CSR_double_row(A,i), k)
#define CSR_double_row_nz(A, i)	 SPVEC_double_nz(CSR_double_row(A,i))
#define CSR_double_nz(A)		((A)->nz_)
#define CSR_double_M(A)		((A)->M_)
#define CSR_double_N(A)		((A)->N_)



CSR_double  *CSR_double_new( int M, int N);

void CSR_double_delete(CSR_double *A);

int CSR_double_insert_entry(CSR_double *A, double val, int i, int j);

int CSR_double_insert_entries(CSR_double *A, int nz, const double *val, 
			const int *I, const int *J);
 
int CSR_double_insert_entries_1(CSR_double *A, int nz, const double *val, 
		const int *I, const int *J);
 
int CSR_double_insert_row(CSR_double *A, int row, int nz, const double *val, 
								const int *index);

int CSR_double_insert_row_1(CSR_double *A, int row, int nz, const double *val, 
								const int *index);

int CSR_double_insert_col(CSR_double *A, int col, int nz, const double *val, 
				const int *index);

int CSR_double_insert_clique(CSR_double *A, int k, int l, const double *val, 
		int row_stride, int col_stride, const int *I, const int *J);

int CSR_double_insert_block(CSR_double *A, const double *val,  int row_stride, 
		int col_stride, int i, int j, int k, int l);

void CSR_double_trim(CSR_double *A);

void CSR_double_end(CSR_double *A);

void CSR_double_print(const CSR_double *A, const char *format);


/*-------------------------*/
/*  Computational Routines */
/*-------------------------*/


void CSR_double_mv(double alpha, const CSR_double *A,  const double *x, 
		int incx, double *y, int incy);

void CSR_double_mm(double alpha, const CSR_double *A,  const double *B, int K,
		int Brow_stride, int Bcol_stride, double *C, 
		int Crow_stride, int Ccol_stride);

void CSR_double_mtv(double alpha, const CSR_double *A,  const double *x, 
		int incx, double *y, int incy);

void CSR_double_mtm(double alpha, const CSR_double *A,  const double *B, 
		int K, int Brow_stride, int Bcol_stride, double *C, 
		int Crow_stride, int Ccol_stride);



#endif
