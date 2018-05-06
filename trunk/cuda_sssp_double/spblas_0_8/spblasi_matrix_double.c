
/**

	BLAS Internal Sparse Matrix Object:

	Holds a generic Sparse BLAS matrix representation, which includes
	a possible variable-block or constant-block format.  Also includes
	a possible separate matrix diagonal for symmetric or triangular
	systems.  Also includes the an encoding of the Sparse BLAS
	matrix properties, as well as an internal state indicator 
	to tell if the matrix has been initalized, is being built, or
	is ready for BLAS computations.


	Example usage:  

<code>
	SPBLASI_Matrix_prop_type prop;

	SPBLASI_set_zero_base(&p);
	SPBLASI_set_real(&p);
	SPBLASI_set_general(&p);
	SPBLASI_set_double_precision(&p);

	SPBLASI_Matrix  A = SPBLASI_Matrix_double_new(M, N, prop);
</code>

	(at this point one could still modify the matrix properties
	further.)

<code>
	SPBLASI_set_real( & SPBLAS_Matrix_prop_type(A) );
	SPBLASI_Matrix_double_prepare_for_insert(A);

	for (i=0; i<nz; i++)
		SPBLASI_Matrix_double_insert_entry(A, Val[i], I[i], J[i]);

	SPBLASI_Matrix_double_end_construction(A);
</code>

	or
<code>
	SPBLASI_Matrix_double A = SPBLASI_Matrix_double_new(M, N);
	SPBLASI_Matrix_set_zero_base(A);
	SPBLASI_Matrix_set_real(A);
	SPBLASI_Matrix_set_general(A);
	SPBLASI_Matrxi_set_double_precision(A);

	SPBLASI_Matrix_double_prepare_for_insert(A);

	...
</code>



*/


#include <stdio.h>
#include <stdlib.h>
#include "blas_malloc.h"
#include "spblasi_matrix.h"
#include "spblasi_matrix_double.h"
#include "spblasi_error.h"
#include "microblas_double.h"

	
int SPBLASI_Matrix_double_insert_entry( SPBLASI_Matrix *A, 
			double val, int i, int j)
{
	SPBLASI_Matrix_prop_type p;

	BLAS_ASSERT_RETURN( A!=NULL, -1);
	BLAS_ASSERT_RETURN( i < SPBLASI_Matrix_M(A), -1);
	BLAS_ASSERT_RETURN( j < SPBLASI_Matrix_N(A), -1);
	BLAS_ASSERT_RETURN( SPBLASI_Matrix_state(A) != unused, -1);
	BLAS_ASSERT_RETURN( SPBLASI_Matrix_state(A) != valid, -1);

	/*
	* if this is the first call to an "insert" routine, then
	* signify that the "setup-phase" of creation is now over,
	* and we are entering the "instertion-phase".
	*/

#if 0
	if (SPBLASI_Matrix_state(A) == initialized)
	{

		fprintf(stderr, "DIAG: entering diagonal initialization. \n");

		if (SPBLASI_is_symmetric(SPBLASI_Matrix_prop(A)) ||
		    SPBLASI_is_upper_triangular(SPBLASI_Matrix_prop(A)) ||
		    SPBLASI_is_lower_triangular(SPBLASI_Matrix_prop(A)) )
		{	
			/* initialize diagonal */

			int i = 0;
			int minMN = 0;
			
			minMN = SPBLASI_Matrix_M(A) < SPBLASI_Matrix_N(A) ? 
			SPBLASI_Matrix_M(A) : SPBLASI_Matrix_N(A)  ;

			SPBLASI_Matrix_diag(A) = (void *)
				blas_malloc(minMN * sizeof(double));

			printf("DIAG: initalized diag...\n");

			for (i=0; i<minMN; i++)
				SPBLASI_Matrix_double_diag(A)[i] = 0.0;


		}    
		SPBLASI_Matrix_state(A) = open;
	}
#endif

	BLAS_ASSERT_RETURN( SPBLASI_Matrix_state(A) == open, -1);

	p = SPBLASI_Matrix_prop(A);
	if (SPBLASI_is_one_base(p))
	{
		i--;
		j--;
	}


	if (SPBLASI_is_symmetric(p))
	{
		if (i==j)
		{
			SPBLASI_Matrix_double_diag(A)[i] += val;
			return 0;
		}
		/* 
	 	* if off-diagonal, store lower triangular portion only.
		*/
		else if (i < j)
			return CSR_double_insert_entry(
					SPBLASI_Matrix_double_CSR(A), val, j, i);
		else if (i > j)
			return CSR_double_insert_entry(
					SPBLASI_Matrix_double_CSR(A), val, i, j);

	}
	else if (SPBLASI_is_lower_triangular(p))
	{
		if (i > j)
			return CSR_double_insert_entry(
					SPBLASI_Matrix_double_CSR(A), val, j, i);
		else if (i==j)
		{
			SPBLASI_Matrix_double_diag(A)[i] += val;
			return 0;
		}
		else
			/* 
			 * entry is not in lower-part of matrix.
			*/
			return -1;
	}
	else if (SPBLASI_is_upper_triangular(SPBLASI_Matrix_prop(A)))
	{
		if (i < j)
			return CSR_double_insert_entry(
					SPBLASI_Matrix_double_CSR(A), val, i, j);
		else if (i==j)
		{
			SPBLASI_Matrix_double_diag(A)[i] += val;
			return 0;
		}
		else
			/* 
			 * entry is not in upper-part of matrix.
			*/
			return -1;

	}
	
	return CSR_double_insert_entry(SPBLASI_Matrix_double_CSR(A), val, i, j);
}


int SPBLASI_Matrix_double_prebuild_diag(SPBLASI_Matrix *S)
{
	/* we can assume that S is a square matrix (M==N) */

	int N = S->M_;  /* or S->N_ */
	int i=0;

    S->diag_ = (double*) blas_malloc(sizeof(double)* (N));
	BLAS_ASSERT_RETURN( SPBLASI_Matrix_double_diag(S) != NULL, -1);
	S->minMN_ = N;
	for (i=0; i<N; i++)
		SPBLASI_Matrix_double_diag(S)[i] = 0.0;

	return 0;
}

/*
 * build_diag() is an internal function, not used elsewhere. 
 * (Not included in SPBLASI_Matrix_double.h header file.)
 *
 * At this point, we can assume that matrix is either
 * symmetric or triangular.
 *
*/

int SPBLASI_Matrix_double_build_diag(SPBLASI_Matrix *S)
{
    int i=0;
	CSR_double *A = NULL;
	int M = 0;
	int N = 0;
    int minMN = 0;

	BLAS_ASSERT_RETURN( S != NULL, -1);
	BLAS_ASSERT_RETURN( SPBLASI_Matrix_state(S) != valid, -1);

	A = SPBLASI_Matrix_double_CSR(S);
	M = SPBLASI_Matrix_M(S);
	N = SPBLASI_Matrix_N(S);
    minMN = M < N ?  M : N;


    S->diag_ = (double*) blas_malloc(sizeof(double)* (minMN));
	S->minMN_ = minMN;

	BLAS_ASSERT_RETURN( SPBLASI_Matrix_double_diag(S) != NULL, -1);

    /* initialize diagonal elements */
    for (i=0; i<minMN; i++)
        SPBLASI_Matrix_double_diag(S)[i] = 0.0;

    /* for each row in A */
    for (i=0; i< A->M_; i++)
    {
        int nz = CSR_double_row_nz(A, i);
        int t;

        /* identify diagonals and sum them in place */
        for (t=0; t < nz; t++)
        {
            /* if element is a diagonal, */
            if (CSR_double_index(A, i, t) == i )
            {
                /* add to previous diagonal, and zero-out current slot */
                SPBLASI_Matrix_double_diag(S)[i] += CSR_double_val(A, i, t);
                CSR_double_val(A, i, t) = 0.0;
            }
        }

    }

    return 0;
}


/*
* this is called before the first insert, after all of the matrix
* properties have been set.  Check for non-square symmetric and
* hermitian matrices.  Check that hermitian matrices are complex.
* If things look OK, set the matrix state to "open" for new insertions.
*/
int SPBLASI_Matrix_double_prepare_for_insert(SPBLASI_Matrix *A)
{

    SPBLASI_Matrix_prop_type p = SPBLASI_Matrix_prop(A);
    int M = SPBLASI_Matrix_M(A);
    int N = SPBLASI_Matrix_N(A);
    int symm = SPBLASI_is_symmetric(p);
    int triang = SPBLASI_is_triangular(p);

    BLAS_ASSERT_RETURN(SPBLASI_Matrix_state(A) == initialized, -1);

    if (symm || triang)
	{
        BLAS_ASSERT_RETURN( M == N, -1);
		BLAS_ASSERT_RETURN(SPBLASI_Matrix_double_prebuild_diag(A) == 0, -1);
	}

    SPBLASI_Matrix_state(A) = open;
    return 0;
}


int SPBLASI_Matrix_double_end_construction(SPBLASI_Matrix *A)
{

//	printf("DIAG: entered end_construction.\n");

	BLAS_ASSERT_RETURN(A!=NULL, -1);

	/* if already closed (previous call to this routine)
	*  perform no-op.
	*/
	if(SPBLASI_Matrix_state(A) == valid)
		return 0;


#if 0
	/* this was from before when we were creating the diagonal
	* after insertions.  But now we create it before the
	* first insertion. (See prepare_for_insertion().
	*/
	if (!SPBLASI_is_general(SPBLASI_Matrix_prop(A)))
	{
		/* if not general, must be symmetric or triangular and
		*  hence, square.  
		*/
		if (SPBLASI_Matrix_M(A) != SPBLASI_Matrix_N(A))
		{
			SPBLASI_Matrix_state(A) =  unavailable;
			return -1;
		}
		BLAS_ASSERT_RETURN(SPBLASI_Matrix_double_build_diag(A) != 0, -1);
	}
#endif	

	CSR_double_trim(SPBLASI_Matrix_double_CSR(A));

	SPBLASI_Matrix_state(A) = valid;

	return 0;
}

int SPBLASI_Matrix_double_delete(SPBLASI_Matrix *A)
{
	BLAS_ASSERT_RETURN( A!=NULL, -1);

	CSR_double_delete(SPBLASI_Matrix_double_CSR(A));
	if (A->diag_ != NULL)
		blas_free(A->diag_);

	if (A->K_ != NULL)
		blas_free(A->K_);

	if (A->L_ != NULL)
		blas_free(A->L_);

	return 0;
}


SPBLASI_Matrix *SPBLASI_Matrix_double_new(int M, int N, 
										SPBLASI_Matrix_prop_type prop)
{

  	SPBLASI_Matrix *A = (SPBLASI_Matrix *) 
				blas_malloc(sizeof(SPBLASI_Matrix));

	BLAS_ASSERT_RETURN(A!=NULL, NULL);

	A->csr_ = CSR_double_new(M,N);

	A->prop_ = prop;
	A->state_ = initialized;	
	SPBLASI_set_real(& SPBLASI_Matrix_prop(A));
	SPBLASI_set_double_precision(& SPBLASI_Matrix_prop(A));

	A->M_ = M;
	A->N_ = N;
	A->k_ = 0;
	A->l_ = 0;
	A->Mb_= 0;
	A->Nb_ = 0;
	A->K_ = NULL;
	A->L_ = NULL;
	A->diag_ = NULL;
	A->minMN_ = M < N ? M : N;

	return A;
}



/** block and variable_block STILL NEED TO BE FIXED! **/

SPBLASI_Matrix *SPBLASI_Matrix_double_new_block(int Mb, int Nb, int k, int l, 
						SPBLASI_Matrix_prop_type prop)
{
	int M = Mb * k;
	int N = Nb * l;

  	SPBLASI_Matrix *A = (SPBLASI_Matrix *) 
				blas_malloc(sizeof(SPBLASI_Matrix));

	BLAS_ASSERT_RETURN(A!=NULL, NULL);

	A->csr_ = CSR_double_new(M,N);

	A->prop_ = prop;
	SPBLASI_set_real(& SPBLASI_Matrix_prop(A));
	SPBLASI_set_double_precision(& SPBLASI_Matrix_prop(A));

	A->state_ = initialized;	
	A->k_ = k;
	A->l_ = l;
	A->Mb_= Mb;
	A->Nb_ = Nb;
	A->K_ = NULL;
	A->L_ = NULL;
	A->diag_ = NULL;
	A->minMN_ = 0;

	return A;
}


SPBLASI_Matrix *SPBLASI_Matrix_double_new_variable_block(int Mb, int Nb, 
		const int *K, const int *L, SPBLASI_Matrix_prop_type prop)
{
	int M = 0;
	int N = 0;
	int i = 0;

  	SPBLASI_Matrix *A = (SPBLASI_Matrix *) 
				blas_malloc(sizeof(SPBLASI_Matrix));

	BLAS_ASSERT_RETURN( A!=NULL, NULL);

	if (Mb == 0 || Nb == 0)
		return NULL;

	for (i=0; i<Mb; i++)
		M += K[i];

	for (i=0; i<Nb; i++)
		N += L[i];

	A->K_ = (int *) blas_malloc((M+1) * sizeof(int));
	BLAS_ASSERT_RETURN(A->K_ != NULL, NULL);

	A->L_ = (int *) blas_malloc((N+1) * sizeof(int));
	BLAS_ASSERT_RETURN(A->L_ != NULL, NULL);



	/* now go through and create a partial sum of block sizes */
	A->K_[0] = 0;
	A->L_[0] = 0;

	for (i=1; i<=Mb; i++)
		A->K_[i] = A->K_[i-1] + K[i-1];

	for (i=1; i<=Nb; i++)
		A->L_[i] = A->L_[i-1] + L[i-1];


	SPBLASI_Matrix_CSR(A) = CSR_double_new(M,N);
	BLAS_ASSERT_RETURN(A->csr_ != NULL, NULL);

	SPBLASI_Matrix_prop(A) = prop;

	SPBLASI_Matrix_state(A) = initialized;	
	SPBLASI_Matrix_k(A) = 0;
	A->l_ = 0;
	A->Mb_= Mb;
	A->Nb_ = Nb;
	A->diag_ = NULL;
	A->minMN_ = 0;

	return A;
}

int SPBLASI_Matrix_double_usmm(enum blas_order_type order, 
			enum blas_trans_type trans, 
			int k, double alpha, const SPBLASI_Matrix *A, 
			const double *B, int lbB, 
			const double *C, int ldC)
{
	return 0;
}



int SPBLASI_Matrix_double_usmv(
			enum blas_trans_type trans, 
			double alpha, const SPBLASI_Matrix *A, 
			const double *x, int incx, 
			double *y, int incy)
{
	BLAS_ASSERT_RETURN(A!=NULL, -1);

	BLAS_ASSERT_RETURN(SPBLASI_Matrix_get_state(A) == valid, -1);

	if (SPBLASI_is_symmetric(A->prop_))
	{
		/*
		 * perform y = alpha(L + L' + D)x + y as
		 *
		 * (1) y = alpha*Lx + y
		 * (2) y = alpha*L'x + y
		 * (3) y = alpha*Dx + y
		 *
		*/

		CSR_double_mv(alpha, SPBLASI_Matrix_double_CSR(A), x, incx, y, incy);
		CSR_double_mtv(alpha, SPBLASI_Matrix_double_CSR(A), x, incx, y, incy);
		ublas_diagmult_double(SPBLASI_Matrix_minMN(A), alpha, 
			SPBLASI_Matrix_diag(A), x, incx, y, incy);
	}
	else if(SPBLASI_is_upper_triangular(A->prop_) ||
			SPBLASI_is_lower_triangular(A->prop_ ))
	{
		if (trans == blas_no_trans)  
			CSR_double_mv(alpha, SPBLASI_Matrix_double_CSR(A), 
						x, incx, y, incy);
		else
			CSR_double_mtv(alpha, SPBLASI_Matrix_double_CSR(A), 
						x, incx, y, incy);

		ublas_diagmult_double(SPBLASI_Matrix_minMN(A), alpha, 
				SPBLASI_Matrix_diag(A), x, incx, y, incy);
	}
	else
	{
		if (trans == blas_no_trans)  
			CSR_double_mv(alpha, SPBLASI_Matrix_double_CSR(A), 
				x, incx, y, incy);
		else
			CSR_double_mtv(alpha, SPBLASI_Matrix_double_CSR(A), 
				x, incx, y, incy);
			
	}
		
	
	return 0;
}

int SPBLASI_Matrix_double_ussm(enum blas_order_type order, 
			enum blas_trans_type trans, int k, double alpha, 
			const SPBLASI_Matrix *A, 
			double *B, int ldB)
{
	return 0;
}
	
/*-----------------------------------------------------------------*/

int SPBLASI_Matrix_double_insert_entries(SPBLASI_Matrix *A, int nz, 		
			const double *Val, const int *I, const int *J)
{
	int i=0;
	for (i=0; i<nz; i++)
		if (SPBLASI_Matrix_double_insert_entry(A, Val[i], I[i], J[i]) < 0)
			return -1;

	return 0;
}






int SPBLASI_Matrix_double_insert_col( SPBLASI_Matrix *A, int col, int nz, 
		const double *Val, const int *I)
{
	if (A==NULL)
		return -1;

	if ( SPBLASI_Matrix_state(A) != open)
	{
		if (SPBLASI_Matrix_state(A) == initialized)
			SPBLASI_Matrix_state(A) = open;
		else
		  return -1;
	}

	if (SPBLASI_Matrix_state(A) != open)
		return -1;

	if (SPBLASI_is_one_base(A->prop_))
	{
		int *Im1 =  (int *) blas_malloc(sizeof(int)*nz);
		int i=0;
		int ret_value = 0;

		for (i=0; i<nz; i++)
		{
			Im1[i] = I[i]-1;
		}
		
		ret_value =  CSR_double_insert_col(SPBLASI_Matrix_double_CSR(A), nz, 
			col-1, Val, Im1);

		blas_free(Im1);

		return ret_value;
	}	
	else
	{
		return CSR_double_insert_col(SPBLASI_Matrix_double_CSR(A), nz, 
			col, Val, I);
	}
}


int SPBLASI_Matrix_double_insert_row( SPBLASI_Matrix *A, int row, int nz, 
		const double *Val, const int *J)
{
	if (A==NULL)
		return -1;

	if ( SPBLASI_Matrix_state(A) != open)
	{
		if (SPBLASI_Matrix_state(A) == initialized)
			SPBLASI_Matrix_state(A) = open;
		else
		  return -1;
	}

	if (A==NULL || SPBLASI_Matrix_state(A) != open)
		return -1;

#if 0
	if (SPBLASI_is_one_base(A->prop_))
	{
		int *Jm1 =  (int *) blas_malloc(sizeof(int)*nz);
		int i=0;
		int ret_value = 0;

		for (i=0; i<nz; i++)
		{
			Jm1[i] = J[i]-1;
		}
		
		ret_value =  CSR_double_insert_row(SPBLASI_Matrix_double_CSR(A), 
				nz, row-1, Val, Jm1);

		blas_free(Jm1);

		return ret_value;
	}	
#endif
	if (SPBLASI_is_one_base(SPBLASI_Matrix_prop(A)))
	{
		return CSR_double_insert_row_1(SPBLASI_Matrix_double_CSR(A), nz, 
			row-1, Val, J);
		
	}
	else
	{
		return CSR_double_insert_row(SPBLASI_Matrix_double_CSR(A), nz, 
			row, Val, J);
	}
}

int SPBLASI_Matrix_double_insert_block(SPBLASI_Matrix *A, const double *Val, 
		int row_stride, int col_stride, int bi, int bj)
{

	if (A==NULL)
		return -1;

	if ( SPBLASI_Matrix_state(A) != open)
	{
		if (SPBLASI_Matrix_state(A) == initialized)
			SPBLASI_Matrix_state(A) = open;
		else
		  return -1;
	}

	if (A==NULL || SPBLASI_Matrix_state(A) != open)
		return -1;


	if ( SPBLASI_is_one_base(SPBLASI_Matrix_prop(A)))
	{
		bi--;
		bj--;
	}

	if ( SPBLASI_Matrix_has_variable_block(A))
	{
		return CSR_double_insert_block(SPBLASI_Matrix_double_CSR(A), 
			Val, row_stride, col_stride, A->K_[bi], A->L_[bj], 
			A->K_[bi+1]-A->K_[bi], A->L_[bj+1]-A->L_[bj] );
	}
	else if (SPBLASI_Matrix_has_block(A))
	{
		return CSR_double_insert_block( SPBLASI_Matrix_double_CSR(A), 
			Val, row_stride, col_stride, A->k_*bi, A->l_*bj, A->k_, A->l_ );
	}
	else
		return -1;

}


/*
* return 1 if any diagonals are zero; 0  otherwise.
*
*/
int	SPBLASI_diag_has_zero(int N, const double *x)
{
	int i;

	for (i=0;i<N; i++)
		if (x[i] == 0.0)
			return 1;

	return 0;
}

/**
* returns 0 upon success, otherwise -1.
*
*/
int SPBLASI_Matrix_double_ussv( enum blas_trans_type trans, 
  			double alpha, 
			const SPBLASI_Matrix *A, 
			double *x, int incx)
{
	int N = SPBLASI_Matrix_N(A);
	int i=0;
	int ii=0;
	int j=0;
	int jj=0;

	const CSR_double *S = NULL;

	if (alpha == 0.0)
		return -1;

	/* check for any zeros on the diagonal */
	BLAS_ASSERT_RETURN(! SPBLASI_diag_has_zero(N, A->diag_), -1);

	N = SPBLASI_Matrix_N(A);
	S = SPBLASI_Matrix_double_CSR(A);

	if (trans == blas_no_trans)
	{
		if (SPBLASI_is_lower_triangular(A->prop_))
		{
			for (i=0, ii=0; i<N; i++, ii += incx)
				x[ii] = (x[ii] - SPVEC_double_dot(CSR_double_row(S,i), x, incx))
									/SPBLASI_Matrix_double_diag(A)[i];
			if (alpha != 1.0)
			{
				double one_over_alpha = 1.0 / alpha;

				for (ii=(N-1)*incx; 0 <= ii; ii -=incx)
					x[ii] *= one_over_alpha;
			}

			return 0;
		}
		else if (SPBLASI_is_upper_triangular(A->prop_))
		{
			for (i=N-1, ii=(N-1)*incx; 0<=i; i--, ii -= incx)
				x[ii] = (x[ii]-SPVEC_double_dot(CSR_double_row(S,i), x, incx))
										/ SPBLASI_Matrix_double_diag(A)[i];
			if (alpha != 1.0)
			{
				double one_over_alpha = 1.0 / alpha;

				for (ii=(N-1)*incx; 0 <= ii; ii -=incx)
					x[ii] *= one_over_alpha;
			}
			
			return 0;
		}
	}
	else if (trans == blas_trans)
	{
		if (SPBLASI_is_lower_triangular(A->prop_))
		{
			for (j=N-1, jj=(N-1)*incx; 0<=j; j--, jj -= incx)
			{
				double tmp = -(x[jj] /= SPBLASI_Matrix_double_diag(A)[j]);
				SPVEC_double_axpy(tmp, CSR_double_row(S,j), x, incx);
				
			}
			if (alpha != 1.0)
			{
				double one_over_alpha = 1.0 / alpha;
				for (jj=(N-1)*incx; 0 <= jj; jj -=incx)
					x[jj] *= one_over_alpha;
			}


			return 0;
		}
		else if (SPBLASI_is_upper_triangular(A->prop_))
		{
			
			for (j=0, jj=0; j<N; j++, jj += incx)
			{
				double tmp = -(x[jj] /= SPBLASI_Matrix_double_diag(A)[j]);
				SPVEC_double_axpy(tmp, CSR_double_row(S,j), x, incx);
			}
			if (alpha != 1.0)
			{
				double one_over_alpha = 1.0 / alpha;
				for (jj=(N-1)*incx; 0 <= jj; jj -=incx)
					x[jj] *= one_over_alpha;
			}
			return 0;
		}

	}


	return -1;
}

/* Diagnostic tools */

void SPBLASI_Matrix_double_print(const SPBLASI_Matrix *A, const char *format)
{

	CSR_double_print(SPBLASI_Matrix_double_CSR(A), format);
	if (SPBLASI_Matrix_double_diag(A) != NULL)
	{
		int M = SPBLASI_Matrix_M(A);
		int N = SPBLASI_Matrix_N(A);
		int minMN = M < N ? M : N;
		int i=0;

		/* print out diagonals */
		printf("Diagonals: \n");
		for (i=0; i<minMN; i++)
		{
			printf("%8d %8d  ", i, i);
			printf( format, SPBLASI_Matrix_double_diag(A)[i]);
			printf("\n");
		}
		printf("\n");
	}
}


