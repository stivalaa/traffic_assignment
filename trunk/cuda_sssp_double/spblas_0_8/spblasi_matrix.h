
#ifndef SPBLAS_MATRIX_H
#define SPBLAS_MATRIX_H

/* #include "blas_enum.h"  */
#include "spblasi_matrix_prop.h"
#include "spblasi_matrix_state.h"
/*
* SPBLASI_Matrix defines the generic combined structure to support
* the various formats of sparse matrices used in this implementation
* of the Sparse BLAS.
*
*/

/*
 * Matrix states:
 *
 * (1) unused: this the unintialized state.
 * (2) initalized: matrix dimensions have been specified; matrix
 *				properties can be set via ussp(); matrix elements
 *				have not been filled in (i.e. no calls to insert
 *				yet.)
 * (3) open:  first insert() is called; properties may no longer
 *				be modified.
 * (4) valid: end (construction) has been called; matrix can now
 *				be used in computational routines.
 *
*/


#define SPBLASI_Matrix_state(A)		((A)->state_) 
#define SPBLASI_Matrix_prop(A)       ((A)->prop_)
#define SPBLASI_Matrix_M(A)			((A)->M_)	
#define SPBLASI_Matrix_N(A)			((A)->N_)	

#define SPBLASI_Matrix_Mb(A)		((A)->Mb_)
#define SPBLASI_Matrix_Nb(A)		((A)->Nb_)
#define SPBLASI_Matrix_K(A)		((A)->K_)
#define SPBLASI_Matrix_L(A)		((A)->L_)
#define SPBLASI_Matrix_k(A)		((A)->k_)
#define SPBLASI_Matrix_l(A)		((A)->l_)
#define SPBLASI_Matrix_has_block(A)             ((A)->Mb_ != 0 )
#define SPBLASI_Matrix_has_variable_block(A)    ((A)->K_ != NULL)

#define SPBLASI_Matrix_diag(A)	((A)->diag_)
#define SPBLASI_Matrix_CSR(A)	((A)->csr_)
#define SPBLASI_Matrix_minMN(A)  ((A)->minMN_)

/* shortcuts for setting SPBLASI_Matrix properties */

#define SPBLASI_Matrix_is_zero_base(A)  \
			SPBLASI_is_zero_base(SPBLASI_Matrix_prop(A))
#define SPBLASI_Matrix_is_one_base(A) \
			SPBLASI_is_one_base(SPBLASI_Matrix_prop(A))

#define SPBLASI_Matrix_set_general(A) \
			SPBLASI_set_general( & SPBLASI_Matrix_prop(A) )
#define	SPBLASI_Matrix_set_complex(A) \
			SPBLASI_set_complex( & SPBLASI_Matrix_prop(A) )
#define SPBLASI_Matrix_set_lower_symmetric(A) \
			SPBLASI_set_lower_symmetric( & SPBLASI_Matrix_prop(A) )

typedef struct
{
	enum  SPBLASI_Matrix_state_type state_;
	SPBLASI_Matrix_prop_type prop_;


	/*
	 * optional block parameters, active only when 
	 * matrix was constructed via block_begin or varible_block_begin.
	 *
	 * If Mb_ and Nb_ are nonzero, then blocking parameters are not used.
	 * If K or L is NULL, then constant block sizes are used.
	*/
	int M_, N_;
	int Mb_, Nb_;
	int k_, l_;
	int *K_, *L_;
	int minMN_;


	/*
	* the actual compressed-row matrix in 4 possible date types:
	* CSR_double, CSR_float, CSR_fcomplex, CSR_dcomplex.
	*
	*/
	void *csr_;

	/* 
	 * The (optional) diagonal of the matrix.  This is used only
	 * when the matrix is symmetric or triangular.  It is of the
	 * the same floating-point type as csr_;
	 *
	*/

	void *diag_;	

	

} SPBLASI_Matrix;



enum SPBLASI_Matrix_state_type SPBLASI_Matrix_get_state( 
		const SPBLASI_Matrix *A);

void SPBLASI_Matrix_set_state( SPBLASI_Matrix *A, 
								 enum SPBLASI_Matrix_state_type s);

SPBLASI_Matrix_prop_type SPBLASI_Matrix_get_property_tag(
									const SPBLASI_Matrix *A);

void SPBLASI_Matrix_set_property_tag(SPBLASI_Matrix *A,
									SPBLASI_Matrix_prop_type prop);
											
SPBLASI_Matrix *SPBLASI_Matrix_new(int M, int N, SPBLASI_Matrix_prop_type prop); 
SPBLASI_Matrix *SPBLASI_Matrix_new_block(int Mb, int Nb, int k, int l, 
						SPBLASI_Matrix_prop_type prop);
SPBLASI_Matrix *SPBLASI_Matrix_new_variable_block(int Mb, int Nb, 
		const int *K, const int *L, SPBLASI_Matrix_prop_type prop);


int SPBLASI_Matrix_delete(SPBLASI_Matrix *A);

int SPBLASI_Matrix_prepare_for_insert(SPBLASI_Matrix *A);


/* Diagnostic tools */
void SPBLASI_Matrix_print_info(const SPBLASI_Matrix *A);

#endif
