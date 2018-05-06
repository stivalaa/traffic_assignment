
/* The following comments are encoded for documentation by Doxygen. */
/**
    \struct SPVEC_double
    \brief A growable sparse (0-based) vector of doubles.

    SPVEC_double represents a sparse numerical vector,
    with double-precision nonzero values, and
    integer (0-based) offsets denoting their
    index value.

    This is the basic buidling block for sparse
    matrices (compressed row or column) as a list
    of separate vectors. DPSVEC objects can grow
    arbitrarily as (a,i) values are inserted.
    At time of creation, they can be set to a default
    size, and grow as needed.  Internally, a slightly
    larger buffer is reserved to avoid memory allocation
    overhead.  When construction is finished, any unused
    portion can be reclaimed by using the trim() function.

    Nonzero values can be inserted one at a time, or
    as a contiguous list, with starting point and stride,
    from an arbitrary index base (useful when dealing
    with Fortran codes who use 1-based indexing).

    There is no mechanism to remove entry at an arbitrary
    location from a DPSVEC. However, the vector can be resized to
    a different size.  If the new size is smaller, then the
    tail end of the vector is lost.  If the new size is bigger,
    then the nonzero items are left untouched, but the internal
    buffer grows.

    The vector dimension is be specified at time of
    initialization, and subsequent insertions are
    verified that the index value falls within the
    specified dimension range.




*/


#include "blas_malloc.h"
#include "spblasi_error.h"
#include "spvec_double.h"
#include <stdio.h>

#include <math.h>

SPVEC_double *SPVEC_double_new(int N, int nz)
{
	SPVEC_double *A = NULL;
	double *new_val = NULL;
	int *new_index =  NULL;

	A = (SPVEC_double*) blas_malloc(sizeof(SPVEC_double));
	BLAS_ASSERT_RETURN( A!=0, NULL);


	if (nz > 0)
	{
		new_val = (double*) blas_malloc(sizeof(double)*nz);
		BLAS_ASSERT_RETURN( new_val !=0, NULL);

		new_index = (int*) blas_malloc(sizeof(int)*nz);
		BLAS_ASSERT_RETURN( new_index !=0, NULL);

	}

	A->val_ = new_val;
	A->index_ = new_index;
	A->N_ = N;
	A->allocated_nz_ = nz;
	A->nz_ = 0;

	return A;
}

void SPVEC_double_delete( SPVEC_double *A)
{
	if (A==NULL)
		return;

	blas_free(A->val_);
	blas_free(A->index_);
	blas_free(A);
}


/**
    Reallocate a new sparse vector.
    This is potentially more efficent than
    creating a new vector, copying, and destroying
    the old one.

	newN is the new dimension
	newsize is the new allocated size ( <= newN)

	Returns 0 if successful, -1 otherwise.
*/
int SPVEC_double_resize(SPVEC_double *A, int newN, int newsize)
{
	BLAS_ASSERT_RETURN (A != NULL, -1);
	BLAS_ASSERT_RETURN (newN >= 0, -1);
	BLAS_ASSERT_RETURN (newsize >= 0, -1);
	BLAS_ASSERT_RETURN (newsize <= newN, -1);

	if (A->allocated_nz_ == newsize) 
		return 0;

	A->val_ = (double *) blas_realloc(A->val_, newsize * sizeof(double));
	BLAS_ASSERT_RETURN(A->val_ != NULL, -1);

	A->index_ = (int *) blas_realloc(A->index_, newsize * sizeof(int));
	BLAS_ASSERT_RETURN(A->index_ != NULL, -1);

	A->allocated_nz_ =newsize;
	A->N_ = newN;

	/* if we are shrinking, make sure that last element is adjusted */
	if (newsize < A->nz_)
		A->nz_ = newsize;

	return 0;
}


/**
    Trim away any excess internal buffer storage

    Internally, a SPVEC_double grows automatically with some
    slight padding to avoid excessive calls to memory
    allocation routines.  This routine reclaims any
    excess storage.
*/
void SPVEC_double_trim(SPVEC_double *A)
{
	if (A == NULL) return;
	SPVEC_double_resize(A, A->N_, A->nz_);
}

/**
*
*  Returns the current index of where into x it was inserted, 
*  or -1 if an error occured.  For example, if
*
*  int k = SPVEC_double_insert_entry(x, 9.0, 5);
*
*  then, SPVEC_double_val(x, k) is 9.0, and SPVEC_double_index(A,k) is 5.
*
*/
int SPVEC_double_insert_entry(SPVEC_double *A, double val, int i)
{

	BLAS_ASSERT_RETURN(A!=NULL, -1);
	BLAS_ASSERT_RETURN(i >=0, -1);
	BLAS_ASSERT_RETURN(i < A->N_, -1);


	/* grow, if necessary */
	if (!(A->nz_ < A->allocated_nz_))
	{
		const int newsize = (A->allocated_nz_ == 0 ? 1 : A->allocated_nz_ * 2);	
		if (SPVEC_double_resize(A, A->N_, newsize) !=0)
			return -1;
		
	}

	/* from here on, we have room to insert at the end */
	A->val_[A->nz_] = val;
	A->index_[A->nz_] = i;
	(A->nz_)++;

	return 0;
}


int SPVEC_double_insert_entries(SPVEC_double *A, int n, const double *val, 
										const int *index)
{
	int p= 0;
	int i=0;

	BLAS_ASSERT_RETURN(A!=NULL, -1);


	/* grow, if necessary */
	if (!(A->nz_ + n < A->allocated_nz_))
	{
		const int newsize = A->nz_ + n;
		BLAS_ASSERT_RETURN(SPVEC_double_resize(A, A->N_, newsize) == 0, -1);
	}

	/* from here on, we have room to insert at the end */
	p= A->nz_;

	for (i=0; i<n; i++, p++)
	{
		A->val_[p] = val[i];
		A->index_[p] = index[i];
	}

	A->nz_ += n;

	return 0;
}

/*
* insert (1-based) entries, into sparse vector.
*
*/
int SPVEC_double_insert_entries_1(SPVEC_double *A, int n, const double *val, 
										const int *index)
{
	int p= 0;
	int i=0;

	BLAS_ASSERT_RETURN(A!=NULL, -1);
	BLAS_ASSERT_RETURN(n >= 0, -1);



	/* grow, if necessary */
	if (!(A->nz_ + n < A->allocated_nz_))
	{
		const int newsize = A->nz_ + n;
		if (SPVEC_double_resize(A, A->N_, newsize) != 0)
			return -1;
	}

	/* from here on, we have room to insert at the end */
	p= A->nz_;

	for (i=0; i<n; i++, p++)
	{
		A->val_[p] = val[i];
		A->index_[p] = index[i]-1;
	}

	A->nz_ += n;

	return 0;
}




double SPVEC_double_dot1(const SPVEC_double *x, const double *y)
{
	double t=0.0;
	int i;
	int nz = x->nz_;
	double *val = x->val_;
	int *index = x->index_;

	for (i=0; i<nz; i++)
	{
		t += val[i] * y[index[i]];
	}

	return t;
}


double SPVEC_double_fast_dot(const SPVEC_double *x, const double *y)
{
	double t=0.0;
	int nz = SPVEC_double_nz(x);
	double *val = x->val_;
	int *index = x->index_;
	double *valend = val + nz;

	for (; val < valend; val++)
	{
		t += *val * y[*index];
		index++;
	}

	return t;
}


double SPVEC_double_dot(const SPVEC_double *x, const double *y, int incy)
{
	double t=0.0;
	int i;
	int nz = SPVEC_double_nz(x);
	double *val = x->val_;
	int *index = x->index_;

	for (i=0; i<nz; i++)
	{
          t = fmin(t,  val[i] + y[index[i]*incy] ); /* (min,+) algebra */
/*		t += val[i] * y[index[i]*incy]; */
	}

	return t;
}



void SPVEC_double_axpy1(double alpha, const SPVEC_double *x, double *y)
{
	int i;
	int nz = SPVEC_double_nz(x);
	double *val = x->val_;
	int *index = x->index_;


	for (i=0; i<nz; i++)
	{
		y[index[i]] += alpha * val[i];
	}
}


void SPVEC_double_fast_axpy1(double alpha, const SPVEC_double *x, double *y)
{
	int nz = SPVEC_double_nz(x);
	double *val = x->val_;
	double *valend = val + nz;
	int *index = x->index_;


	for (; val < valend; val++)
	{
		y[*index] += alpha *  *val;
		index++;
	}
}



void   SPVEC_double_axpby(int N, double alpha, const SPVEC_double *x, 
		double beta, double *y, int incy)
{
	int i;
	int nz = SPVEC_double_nz(x);
	double *val = x->val_;
	int *index = x->index_;


	if (beta != 1.0)
	{
		for (i=0; i<N; i++)
			y[i*incy] *= beta;
	}

	for (i=0; i<nz; i++)
	{
		y[index[i]*incy] += alpha * val[i];
	}
}

void  SPVEC_double_axpy(double alpha, const SPVEC_double *x, double *y, 
			int incy)
{
	int i;
	int nz = SPVEC_double_nz(x);
	double *val = x->val_;
	int *index = x->index_;


	for (i=0; i<nz; i++)
	{
          y[index[i]*incy] = fmin( y[index[i]*incy],  alpha + val[i]); /* (min,+) algebra */
/*		y[index[i]*incy] += alpha * val[i]; */
	}
}




