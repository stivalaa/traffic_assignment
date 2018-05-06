
#ifndef SPVEC_double_H
#define SPVEC_double_H 1

/*-------------------------------*/
/*  A "growable" sparse vector   */
/*-------------------------------*/


typedef struct
{
	int N_;
    int allocated_nz_;	/* original value for malloc */
	int nz_;
    double *val_;      /* scalar values */
	int *index_;		/* offsets */
} SPVEC_double;

SPVEC_double *SPVEC_double_new(int N, int nz);
void  SPVEC_double_delete(SPVEC_double *A);
int   SPVEC_double_resize(SPVEC_double *A, int newN, int newsize);
void  SPVEC_double_trim(SPVEC_double *A);


int SPVEC_double_insert_entry(SPVEC_double *A, double val, int i);
int SPVEC_double_insert_entries(SPVEC_double *A, int n, const double *val, 
			const int *index);
int SPVEC_double_insert_entries_1(SPVEC_double *A, int n, const double *val, 
			const int *index);

#define SPVEC_double_val(A,i)  		((A)->val_[i])
#define SPVEC_double_index(A,i)		((A)->index_[i])
#define SPVEC_double_nz(A)			((A)->nz_)
#define SPVEC_double_dim(A)			((A)->N_)
#define SPVEC_double_allocated_nz(A)	((A)->allocated_nz_)




/* -------------------------- SPARSE DOT PRODUCT ------------------------ */
/*                                                                        */
/*                              r = dot(x,y);                             */
/*                                                                        */
/* ---------------------------------------------------------------------- */


double SPVEC_double_dot1(const SPVEC_double *x, const double *y);
double SPVEC_double_fast_dot1(const SPVEC_double *x, const double *y);
double SPVEC_double_dot(const SPVEC_double *x, const double *y, int yinc);
void   SPVEC_double_axpy1(double alpha, const SPVEC_double *x, double *y) ;
void   SPVEC_double_fast_axpy1(double alpha, const SPVEC_double *x, double *y) ;
void   SPVEC_double_axpy(double alpha, const SPVEC_double *x, double *y, int incy);
void   SPVEC_double_axpby(int N, double alpha, const SPVEC_double *x, double beta, 
									double *y, int incy);


#endif
