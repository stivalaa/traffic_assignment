#include "blas.h"
#include "spblasi_error.h"
#include "blas_malloc.h"
#include "csr_double.h"

#include <math.h>

/*-----------------------------------------------------------------*/
/*  Compressed Row Sparse Matrix: based on list of sparse vectors  */
/*-----------------------------------------------------------------*/


CSR_double *CSR_double_new( int M, int N)
{
	CSR_double *A = NULL;
	int i =0;

	BLAS_ASSERT_RETURN(M > 0, NULL);
	BLAS_ASSERT_RETURN(N > 0, NULL);

	A = (CSR_double*) blas_malloc(sizeof(CSR_double));
	BLAS_ASSERT_RETURN( A != NULL, NULL);

	A->row_ = (SPVEC_double**) blas_malloc(sizeof(SPVEC_double*) * M);
	BLAS_ASSERT_RETURN( A->row_ != 0, NULL);


	/* now, allocate each sparse vector */
	for (i=0; i<M; i++)
	{
		A->row_[i] = SPVEC_double_new(N, 0);
	}


	A->M_ = M;
	A->N_ = N;
	A->nz_ = 0;

	return A;
}

/**
*
*   insert a new nonzero entry into location A(i,j).
*   
*
*	@return 0 if sucessful, -1 otherwise.
*/
int CSR_double_insert_entry(CSR_double *A, double val, int i, int j)
{


	BLAS_ASSERT_RETURN( A!= NULL, -1);
	BLAS_ASSERT_RETURN( i >= 0, -1);
	BLAS_ASSERT_RETURN( i < CSR_double_M(A), -1);
	BLAS_ASSERT_RETURN( j >= 0, -1);
	BLAS_ASSERT_RETURN( j < CSR_double_N(A), -1);


	if (SPVEC_double_insert_entry(A->row_[i], val, j) != 0) 
			return -1;

	A->nz_ ++;
	
	return 0;
}


int CSR_double_insert_row(CSR_double *A, int row, int nz, const double *val, 
						const int *index)
{

 	if (A==NULL || nz < 0 )
		return -1;
 	if ( row < 0 || row >= A->M_ ) 
		return -1;


	if (SPVEC_double_insert_entries(A->row_[row], nz, val, index) != 0)
		return -1;

	A->nz_ += nz;
	return 0;
}

int CSR_double_insert_row_1(CSR_double *A, int row, int nz, const double *val, 
						const int *index)
{

 	if (A==NULL || nz < 0 )
		return -1;
 	if ( row < 0 || row >= A->M_ ) 
		return -1;


	if (SPVEC_double_insert_entries_1(A->row_[row], nz, val, index) != 0)
		return -1;

	A->nz_ += nz;
	return 0;
}


/**
*
*	Insert a list of coordinate entries
*
*	@return -1 if unsucessful. Otherwise, the offset to new location
*
*/
int CSR_double_insert_entries(CSR_double *A, int nz, const double *val, const int *I, 
			const int *J)
{
	int i=0;

 	if (nz < 0)
		return -1;

	for (i=0; i<nz; i++)
	{
		if (CSR_double_insert_entry(A, val[i], I[i], J[i]) !=0)
			return -1;
	}

	return 0;
}

/**
*
*	Insert a list of (1-based) coordinate entries
*
*	@return -1 if unsucessful. Otherwise, the offset to new location
*
*/
int CSR_double_insert_entries_1(CSR_double *A, int nz, const double *val, const int *I, 
			const int *J)
{
	int i=0;

 	if (nz < 0)
		return -1;

	for (i=0; i<nz; i++)
	{
		if (CSR_double_insert_entry(A, val[i], I[i]-1, J[i]-1) !=0)
			return -1;
	}

	return 0;
}


int CSR_double_insert_col(CSR_double *A, int col, int nz, const double *val, 
							const int *index)
{
	int i = 0;

 	if (nz < 0)
		return -1;

	if (col < 0 || col >= A->N_)
		return -1;

	for (i=0; i<nz; i++)
	{
		if (CSR_double_insert_entry(A, val[i], index[i], col) != 0)
			return -1;
	}
	
	return 0;
}


int CSR_double_insert_clique(CSR_double *A, int k, int l, 
		const double *val, int row_stride, int col_stride,
		const int *I, const int *J) 
{
	int i=0, j=0;

	for (i=0; i<k; i++)
	{
		int row = I[i];
		for (j=0; j<l; j++)
		{
		   if (CSR_double_insert_entry(A, val[i*row_stride+j*col_stride], row, J[j]))
			return -1;

		}
	}
	return 0;
}

int CSR_double_insert_block(CSR_double *A, const double *val, int row_stride, 
			int col_stride, int i, int j, int k, int l)
{
	int ii=0, jj=0;

	for (ii=0; ii<k; ii++)
		for (jj=0; jj<l; jj++)
		{
			if (CSR_double_insert_entry(A, val[ii*row_stride + jj*col_stride], 
							i+ii, j+jj))
				return -1;
		}

	return 0;
}


/**
*
* This routine is called when construction is completed, i.e. 
* all insertions are complete.  It basically trims excess storage
* no longer needed (i.e. make the capacity the same as the actual
* size).  
*
*/
void CSR_double_trim(CSR_double *A)
{
	int i=0;

	/* for each row in A */
	for (i=0; i< A->M_; i++)
	{
		SPVEC_double_trim( A->row_[i]);
	}
}


/**
* CSR_double_end() is the same as CSR_double_trim.
*
*/
void CSR_double_end(CSR_double *A)
{
	int i=0;

	/* for each row in A */
	for (i=0; i< A->M_; i++)
	{
		SPVEC_double_trim( A->row_[i]);
	}
}





void CSR_double_delete( CSR_double *A)
{
	int i=0;

	if (A==NULL)
		return ;

	for (i=0; i<A->M_; i++)
	{
		SPVEC_double_delete(A->row_[i]);
	}
	blas_free(A->row_);
	blas_free(A);
}




/*
 * -------------------  SPARSE MATRIX MULTIPLY  ------------------------
 *                  (compressed-row sparse storage)
 *
 *                      y = alpha * A*x + y
 *                      C = alpha * A*B + C
 *
 *  A is sparse (compressed-row), B, and C are dense arrays; x and y
 *  are dense vectors.
 *
 * ----------------------------------------------------------------------
*/

void CSR_double_mv(double alpha, const CSR_double *A,  const double *x, 
	int incx, double *y, int incy)
{
	int M = CSR_double_M(A);
	int i;

    for (i=0; i<M; i++)
      y[i*incy] = fmin( y[i*incy],  alpha + SPVEC_double_dot(CSR_double_row(A,i), x, incx) ); /* (min,+) algebra */
/*         y[i*incy] += alpha * SPVEC_double_dot(CSR_double_row(A,i), x, incx); */
}

/*
 * matrix transpose multiply.
 *
*/
void CSR_double_mtv(double alpha, const CSR_double *A,  const double *x, 
	int incx, double *y, int incy)
{
	int i;
	int M = CSR_double_M(A);

	for (i=0; i<M; i++)
          SPVEC_double_axpy(alpha + x[i], CSR_double_row(A, i), y, incy); /* (min,+) algebra */
/*		SPVEC_double_axpy(alpha * x[i], CSR_double_row(A, i), y, incy); */

}

/* A is MxN, B is NxK, C is MxK */

void CSR_double_mm(double alpha, const CSR_double *A,  const double *B, int K,
		int B_row_stride, int B_col_stride, 
		double *C, int C_row_stride, int C_col_stride)
{
    /* perform K matvecs: ith-column of C gets A*(i-th column of B) */

    int i;
    for (i=0; i<K; i++)
        CSR_double_mv(alpha, A, &B[i*B_row_stride], B_col_stride, 
            &C[i*C_row_stride], C_col_stride);
}

/*
 * Tranpose multiply: A is NxM, B is NxK, C is MxK.
 *
*/
void CSR_double_mtm(double alpha, const CSR_double *A,  const double *B, int K,
		int B_row_stride, int B_col_stride, double *C, 
		int C_row_stride, int C_col_stride)
{
    int i;
    for (i=0; i<K; i++)
        CSR_double_mtv(alpha, A, &B[i*B_row_stride], B_col_stride, 
            &C[i*C_row_stride], C_col_stride);


}

/*
*  I/O diagnostic.
*/

#include <stdio.h>
void CSR_double_print(const CSR_double *A, const char *format)
{
	int i=0;
	int j=0;

	if (format == NULL)
		format = "% 20.18e";

	printf("A is of size (%d x %d) with %d nonzeros\n", CSR_double_M(A),
		CSR_double_N(A), CSR_double_nz(A) );



	for (i=0; i< CSR_double_M(A); i++)
	{
		SPVEC_double *rowi = CSR_double_row(A,i);
		int row_nz = SPVEC_double_nz(rowi);

		/* printf("row[%d], nz=%d: \n",i, row_nz); */

		for (j=0; j< row_nz; j++)
		{
			printf("%8d %8d  ", i, SPVEC_double_index(rowi, j));
			printf( format,  SPVEC_double_val(rowi, j));
			printf("\n");
		}
	}

}


