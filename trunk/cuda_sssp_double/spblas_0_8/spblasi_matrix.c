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


	SPBLASI_Matrix_prop_type prop;

	SPBLASI_set_zero_base(&p);
	SPBLASI_set_real(&p);
	SPBLASI_set_general(&p);
	SPBLASI_set_double_precision(&p);

	SPBLASI_Matrix  A = SPBLASI_Matrix_double_new(M, N, prop);

	(at this point one could still modify the matrix properties
	further.)

	SPBLASI_set_real( & SPBLAS_Matrix_prop_type(A) );
	SPBLASI_Matrix_double_prepare_for_insert(A);

	for (i=0; i<nz; i++)
		SPBLASI_Matrix_double_insert_entry(A, Val[i], I[i], J[i]);

	SPBLASI_Matrix_double_end_construction(A);


	or

	SPBLASI_Matrix_double A = SPBLASI_Matrix_double_new(M, N);
	SPBLASI_Matrix_set_zero_base(A);
	SPBLASI_Matrix_set_real(A);
	SPBLASI_Matrix_set_general(A);
	SPBLASI_Matrxi_set_double_precision(A);

	SPBLASI_Matrix_double_prepare_for_insert(A);

	...




*/


#include <stdio.h>
#include <stdlib.h>
#include "blas_malloc.h"
#include "spblasi_matrix.h"
#include "spblasi_matrix_double.h"
#include "spblasi_error.h"


/* These functions are genreic to any BLASI_Matrix */

enum SPBLASI_Matrix_state_type SPBLASI_Matrix_get_state( const SPBLASI_Matrix *A)
{
	return A->state_;
}

void SPBLASI_Matrix_set_state( SPBLASI_Matrix *A, 
								 enum SPBLASI_Matrix_state_type s)
{
	A->state_ = s;
}

SPBLASI_Matrix_prop_type SPBLASI_Matrix_get_property_tag(
									const SPBLASI_Matrix *A)
{
	return A->prop_;
}

void SPBLASI_Matrix_set_property_tag(SPBLASI_Matrix *A,
									SPBLASI_Matrix_prop_type prop)
{
	A->prop_ = prop;
}



int SPBLASI_Matrix_corrupt(SPBLASI_Matrix *A)
{
	SPBLASI_Matrix_state(A) = unavailable;	
	return -1;
}
/*
* this is called before the first insert, after all of the matrix
* properties have been set.  Check for non-square symmetric and
* hermitian matrices.  Check that hermitian matrices are complex.
* If things look OK, set the matrix state to "open" for new insertions.
*/
int SPBLASI_Matrix_prepare_for_insert(SPBLASI_Matrix *A)
{

	SPBLASI_Matrix_prop_type p = SPBLASI_Matrix_prop(A);
	int M = SPBLASI_Matrix_M(A);
	int N = SPBLASI_Matrix_N(A);
	int symm = SPBLASI_is_symmetric(p);
	int triang = SPBLASI_is_triangular(p);
	int herm = SPBLASI_is_hermitian(p);

	BLAS_ASSERT_RETURN(SPBLASI_Matrix_state(A) == initialized, -1);

	if (symm || triang)
		BLAS_ASSERT_RETURN( M == N, -1);

	if (herm)
		BLAS_ASSERT_RETURN(SPBLASI_is_complex(p), -1);


	SPBLASI_Matrix_state(A) = open;
	return 0;
}



int SPBLASI_Matrix_delete(SPBLASI_Matrix *A)
{
	SPBLASI_Matrix_prop_type p = SPBLASI_Matrix_prop(A);

	if (SPBLASI_is_real(p) && SPBLASI_is_double_precision(p))
			return SPBLASI_Matrix_double_delete(A);
	else
		return -1;
}
	

void SPBLASI_Matrix_print_info(const SPBLASI_Matrix *A)
{

	if (A == NULL)
	{
		printf("[Sparse BLAS matrix is NULL.] ");
		return;
	}

	switch(SPBLASI_Matrix_state(A))
	{
		case unused:		printf("[unused] ");
							return;

		case initialized:	printf("[initalized] ");
							break;

		case open:			printf("[open] ");
							break;

		case valid:			printf("[valid] ");
							break;

		default:			printf("[invalid state] ");
							return;
	}

	printf(" M: %d  N: %d   ", SPBLASI_Matrix_M(A), SPBLASI_Matrix_N(A));

	SPBLASI_Matrix_prop_print( & SPBLASI_Matrix_prop(A) );


}
