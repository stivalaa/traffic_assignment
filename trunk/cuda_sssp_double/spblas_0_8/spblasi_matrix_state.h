#ifndef SPBLASI_MATRIX_STATE_H
#define SPBLASI_MATRIX_STATE_H

/*
 * Matrix states:
 *
 * (1) unused: this the unintialized state.
 *
 * (2) initalized: matrix dimensions have been specified; matrix
 *				properties can be set via ussp(); matrix elements
 *				have not been filled in (i.e. no calls to insert
 *				yet.)
 *
 * (3) open:  first insert() is called; properties may no longer
 *				be modified.
 *
 * (4) valid: end (construction) has been called; matrix can now
 *				be used in computational routines.
 *
 * (5) unavailable   : matrix is in an error state, e.g. when creating
 *              a non-square matrix and then try to make it symmetric
 *              via blas_ussp(A, blas_symmetric);
 *
*/

	
enum SPBLASI_Matrix_state_type {unused, initialized, open, valid, unavailable};


#endif
