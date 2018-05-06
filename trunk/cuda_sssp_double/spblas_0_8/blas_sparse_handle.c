#include <stdlib.h>
#include "blas_sparse.h"
#include "spblasi_error.h"  
#include "spblasi_matrix.h"
#include "spblasi_table.h"


/* returns "1" if matrix A has specified property, 0 otherwise. */

int BLAS_usgp(blas_sparse_matrix A, int pname)
{
	SPBLASI_Matrix *S = SPBLASI_table_get(A);

	BLAS_ASSERT_RETURN( S!=NULL, 0);

								/* Note: these are returns, in the case */
								/* statement, so "break" is not needed.  */
	switch (pname)
	{
		case (blas_complex) :	return SPBLASI_is_complex(S->prop_);
		case (blas_real)	:	return SPBLASI_is_real(S->prop_);
		case (blas_double_precision) :
								return SPBLASI_is_double_precision(S->prop_);

		case (blas_symmetric) :	return SPBLASI_is_symmetric(S->prop_);
		case (blas_lower_triangular) : return 
							SPBLASI_is_lower_triangular(S->prop_);
		case (blas_upper_triangular) : return 
							SPBLASI_is_upper_triangular(S->prop_);
		/* ... */

		case (blas_valid_handle) : return 
							SPBLASI_Matrix_state(S) == valid ?  1 : 0;

		default:				return 0;
	}

}	


int BLAS_ussp(blas_sparse_matrix h, int pname)
{
	SPBLASI_Matrix *S = SPBLASI_table_get(h);

	BLAS_ASSERT_RETURN(S != NULL, -1);

	BLAS_ASSERT_RETURN (SPBLASI_Matrix_state(S) == initialized, -1);
								
								/* Note: these set() routines are macros, 	*/
								/* so do not use a temp variable as an		*/
								/* argument, ie. do use S->prop_ explicitly.*/
	switch (pname)
	{
		case (blas_complex) :	SPBLASI_set_complex( &(S->prop_)); break;
		case (blas_real)	:	SPBLASI_set_real(&(S->prop_)); break;
		case (blas_double_precision) :
								SPBLASI_set_double_precision(&(S->prop_)); 
								break;

		case (blas_symmetric)       :
		case (blas_lower_symmetric) : SPBLASI_set_lower_symmetric(&(S->prop_));
										break;
		case (blas_upper_symmetric) : SPBLASI_set_upper_symmetric(&(S->prop_));
										break;
		case (blas_lower_triangular) : 
								SPBLASI_set_lower_triangular(&(S->prop_));
								break;
		case (blas_upper_triangular) : SPBLASI_set_upper_triangular
			(&(S->prop_));
										break;
		/* ... */

		default:				BLAS_ASSERT_RETURN(0, -1);
	}

	return 0;
}


/**
*  Destroys interal sparse matrix representation pointed to by
*  handle "h".  Tries to reclaim any used space used by the creation
*  and usage of this matrix, although this is not reqired by the
*  standard.  If "h" does not refer to a valid matrix, function
*  returns -1.  Otherwise, it returns 0 if successful.
*
*/
int BLAS_usds(blas_sparse_matrix h)
{
	SPBLASI_Matrix *A = SPBLASI_table_remove(h);

	if (A==NULL) return -1;

	return SPBLASI_Matrix_delete(A);
}







