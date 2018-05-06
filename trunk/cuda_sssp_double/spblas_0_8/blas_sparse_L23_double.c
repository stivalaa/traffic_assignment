#include <stdlib.h>
#include "blas_malloc.h"
#include "blas_sparse.h"
#include "spblasi_error.h"  
#include "spblasi_matrix_prop.h"
#include "spblasi_matrix_double.h"
#include "spblasi_table.h"



/*
* Level 2 & 3 routines
*
*/




static const SPBLASI_Matrix_prop_type C_default_prop =
{
    blasi_zero_base,
    blasi_rowmajor,
    blasi_no_repeated_indices,
    blasi_general,
    blasi_real,
	blasi_double_precision,
    blasi_irregular
};


static const SPBLASI_Matrix_prop_type F77_default_prop =
{
    blasi_one_base,
    blasi_colmajor,
    blasi_no_repeated_indices,
    blasi_general,
    blasi_real,
	blasi_double_precision,
    blasi_irregular
};


/**
 Begin constructing a new M x N sparse matrix. The default 
 characteristics are 1) an index base of 0, 2) row-major ordering
 for dense blocks, 3) allows repeated indices, and 4) has no
 special symmetry properites (i.e. not symmetric).  Any of these
 may be overriden with calls to BLAS_ussp().

	@param M the row dimension of the sparse matrix
	@param N the column dimension of the sparse matrix.
	@return a handle to internal matrix. (An integer offset into
	the global table.)

*/
int  BLAS_duscr_begin(int M, int N )
{

	SPBLASI_Matrix *A = SPBLASI_Matrix_double_new(M, N, C_default_prop);
	return (A==NULL) ? -1 : SPBLASI_table_insert(A);

}

int  BLAS_duscr_block_begin(int Mb, int Nb, int k, int l)
{
	SPBLASI_Matrix *A = SPBLASI_Matrix_double_new_block(Mb, Nb, k, l,
							C_default_prop);
	return (A==NULL) ? -1 : SPBLASI_table_insert(A);
}


int  BLAS_duscr_variable_block_begin(int Mb, int Nb, const int *k, const int *l)
{
	SPBLASI_Matrix *A = SPBLASI_Matrix_double_new_variable_block(Mb, Nb, k, l,
							C_default_prop);

	return (A==NULL) ? -1 : SPBLASI_table_insert(A);
}


/* ----- Fortran 77 entry points ----------- */

int  f_duscr_begin(int M, int N )
{
	SPBLASI_Matrix *A = 
			SPBLASI_Matrix_double_new(M, N, F77_default_prop);

	return (A==NULL) ? -1 : SPBLASI_table_insert(A);

}

int  f_duscr_block_begin(int Mb, int Nb, int k, int l)
{
	SPBLASI_Matrix *A = SPBLASI_Matrix_double_new_block(Mb, Nb, k, l,
							F77_default_prop);

	return (A==NULL) ? -1 : SPBLASI_table_insert(A);
}


int  f_duscr_variable_block_begin(int Mb, int Nb, const int *k, const int *l)
{
	SPBLASI_Matrix *A = SPBLASI_Matrix_double_new_variable_block(Mb, Nb, k, l,
							F77_default_prop);

	return (A==NULL) ? -1 : SPBLASI_table_insert(A);
}

/* ------------------------------------------- */



int BLAS_duscr_insert_entry(blas_sparse_matrix h, double val, int i, int j)
{
	SPBLASI_Matrix *A = SPBLASI_table_get(h);


	BLAS_ASSERT_RETURN( A!=NULL, -1);
	if (SPBLASI_Matrix_state(A) == initialized)
		BLAS_ASSERT_RETURN(SPBLASI_Matrix_double_prepare_for_insert(A)==0, -1);
	BLAS_ASSERT_RETURN (SPBLASI_Matrix_state(A) == open, -1);

	return SPBLASI_Matrix_double_insert_entry(A, val, i, j);

}

int BLAS_duscr_end(int h)
{
	SPBLASI_Matrix *A = SPBLASI_table_get(h);

	BLAS_ASSERT_RETURN( A!=NULL, -1);
	BLAS_ASSERT_RETURN(SPBLASI_Matrix_get_state(A) == open, -1);

	if (SPBLASI_Matrix_double_end_construction(A) != 0)
		return -1;

	SPBLASI_Matrix_state(A) = valid;

	return 0;
}

int BLAS_duscr_insert_entries(blas_sparse_matrix h, int nz, const double *Val,
					const int *I, const int *J)
{
	SPBLASI_Matrix *A = SPBLASI_table_get(h);

	BLAS_ASSERT_RETURN( A!=NULL, -1);
	if (SPBLASI_Matrix_state(A) == initialized)
		BLAS_ASSERT_RETURN(SPBLASI_Matrix_double_prepare_for_insert(A)==0, -1);
	BLAS_ASSERT_RETURN (SPBLASI_Matrix_state(A) == open, -1);

	
	return SPBLASI_Matrix_double_insert_entries(A, nz, Val, I, J);

}


/*-----------------------------------------------------------------*/
#if 0

int BLAS_duscr_insert_col(blas_sparse_matrix h, int column, int nz, const double *Val, 
						const int *I)
{
	SPBLASI_Matrix *A = SPBLASI_table_get(h);

	if (A == NULL || SPBLASI_Matrix_get_state(A) != open)
		return -1;

	return SPBLASI_Matrix_double_insert_col(A, column, nz, Val, I);

}

int BLAS_duscr_insert_row(blas_sparse_matrix h, int row ,   int nz, const double *Val, 
						const int *I)
{
	SPBLASI_Matrix *A = SPBLASI_table_get(h);

	if (A == NULL || SPBLASI_Matrix_get_state(A) != open)
		return -1;

	return SPBLASI_Matrix_double_insert_row(A, row, nz, Val, I);

}


int BLAS_duscr_insert_clique(blas_sparse_matrix h, int k, int l, 
	const double *val, int row_stride, int col_stride, 
	const int *indx, const int *jndx);

int BLAS_duscr_insert_block(blas_sparse_matrix A, const double *val, 
	int row_stride, int col_stride, int bi, int bj);


#endif








/**
*
*   Sparse matrix times dense vector.  This routine computes
*   y = alpha * op(A) * x + beta * y, where op(A) is either A,
*   the transpose of A, and x and y are dense vectors.
*
*   The routine returns 0 if sucessful, -1 otherwise. If A
*   is an invalid handle, the routine returns -1 and does not
*   modify any of its parameters.
*
*
*   @param trans    compute operation with the transpose of
*   A, if trans is set to blas_transpose.  Otherwise, it computes
*   with A untransposed.   (Note: trans set to blas_conj_trans has
*   the same effect as blas_trans for real valued matrices.)
*
*   @param alpha scalar multiplier of A
*   @param A     (input) sparse matrix handle
*   @param x     (input) dense vector,
*   @param incx  increment between elements of x
*   @param beta  scalar multiplier of y
*   @param y     (input/output) dense vector, one output contains
*                   alpha * op(A) * x + beta * y.
*   @param incy  increment (stride) between elemens of y
*
*/
int BLAS_dusmv(enum blas_trans_type transA, double alpha, 
     blas_sparse_matrix A, const double *x, int incx, double *y, int incy)
{
	SPBLASI_Matrix *S = SPBLASI_table_get(A);

	
	BLAS_ASSERT_RETURN( S!=NULL, -1);
	BLAS_ASSERT_RETURN( transA == blas_no_trans || transA == blas_trans, -1 );
	BLAS_ASSERT_RETURN( BLAS_usgp(A, blas_valid_handle) == 1 , -1);
	BLAS_ASSERT_RETURN( x != NULL, -1 );
	BLAS_ASSERT_RETURN( incx != 0, -1);
	BLAS_ASSERT_RETURN( y != NULL, -1);
	BLAS_ASSERT_RETURN( incy != 0, -1);


	return SPBLASI_Matrix_double_usmv(transA, alpha, SPBLASI_table_get(A), 
					x, incx, y, incy);
}

int BLAS_dusmm(enum blas_order_type order, enum blas_trans_type transA, 
	int k, double alpha, blas_sparse_matrix A, const double *B, int ldB,
	double *C, int ldC)
{
	SPBLASI_Matrix *S = SPBLASI_table_get(A);


	return SPBLASI_Matrix_double_usmm(order, transA, k, alpha, S, B, ldB,
				C, ldC);

}




/**
*
*   Sparse matrix triangular solve.  Given a sparse triangular
*   matrix A and scalar alpha, this routine computes
*   x = alpha * op(A)^{-1}x.
*
*   The routine returns 0 if sucessful, -1 otherwise.
*   The matrix handle A must be refer to an upper or lower triangular
*   matrix, otherwise the routine returns -1,  and does not
*   modify any of its parameters.
*
*   @param trans    computes operation with the transpose of
*   A, if trans is set to blas_transpose.  Otherwise, it computes
*   with A untransposed.   (Note: trans set to blas_conj_trans has
*   the same effect as blas_trans for real valued matrices.)
* 
*   @param alpha scalar multiplier of A.  
*
*   @param T     (input) handle representing a triangular sparse matrix.  
*       Note that T must have been created with the blas_triangular property 
*       set, as well as either the blas_upper or blas_lower property to 
*       denote upper or lower triangular structure.
*
*   @param x     (input/output) on input, it contains the original
*                   right-hand-side of linear equation to be solved.
*               On output, it contains the solution.
*   @param incx  stride between entries of x
*
*/
int BLAS_dussv(enum blas_trans_type transA, double alpha, 
     blas_sparse_matrix A, double *x, int incx)
{

	BLAS_ASSERT_RETURN( transA == blas_trans || transA == blas_no_trans, -1 );
	BLAS_ASSERT_RETURN( BLAS_usgp(A, blas_valid_handle), -1 );
	BLAS_ASSERT_RETURN( BLAS_usgp(A, blas_lower_triangular) ||
				 BLAS_usgp(A, blas_upper_triangular), -1);
	BLAS_ASSERT_RETURN( x != NULL, -1 );
	BLAS_ASSERT_RETURN( incx != 0, -1 );
	

	return SPBLASI_Matrix_double_ussv(transA, alpha, SPBLASI_table_get(A), 
					x, incx);


}

/**

   Sparse matrix triangular solve.  Given a sparse triangular
   matrix A, dense rectangular matrix B and scalar alpha, this routine 
   computes
   B = alpha * op(T)^{-1}B.

   The routine returns 0 if sucessful, -1 otherwise.
   The matrix handle A must be refer to an upper or lower triangular
   matrix, otherwise the routine returns -1,  and does not
   modify any of its parameters.

   @param trans    computes operation with the transpose of
   A, if trans is set to blas_transpose.  Otherwise, it computes
   with A untransposed.   (Note: trans set to blas_conj_trans has
   the same effect as blas_trans for real valued matrices.)
 
   @param alpha scalar multiplier of A.  

   @param A     (input) handle representing a triangular sparse matrix.  
       Note that A must have been created with the blas_triangular property 
       set, as well as either the blas_upper or blas_lower property to 
       denote upper or lower triangular structure.

   @param B     (input/output) on input, it contains the original
                   right-hand-side of linear equation to be solved.
               On output, it contains the solution.
   @param ldB	stride between outer indices of B. (i.e. between rows, if 
   					B is column ordered, or between colums if 
					B is row-ordered.)

*/

int BLAS_dussm(enum blas_order_type order, enum blas_trans_type transA,
	int k, double alpha, blas_sparse_matrix A, double *B, int ldB)
{
	SPBLASI_Matrix *S = SPBLASI_table_get(A);

	return SPBLASI_Matrix_double_ussm(order, transA, k, alpha, S, B, ldB);

}

