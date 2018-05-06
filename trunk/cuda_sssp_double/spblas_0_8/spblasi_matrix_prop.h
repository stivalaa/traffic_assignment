
#ifndef SPBLASI_MATRIX_PROP_H_
#define SPBLASI_MATRIX_PROP_H_ 1





typedef struct
{
	unsigned int index_base : 1;		/* 0 if 0-base, 1 for 1-base */
	unsigned int order : 1;				/* 0 for row, 1 for column */
	unsigned int repeated_indices : 1;	/* 0 for false, 1 for true */


	unsigned int structure : 3;			/* 0 = blasi_general */
										/* 1 = blasi_lower_symmetric */
										/* 2 = blasi_upper_symmetric */
										/* 3 = blasi_lower_triangular */
										/* 4 = blasi_upper_triangualr */
										/* 5 = blasi_lower_hermitian */
										/* 6 = blasi_upper_hermitian */

	unsigned int scalar_type : 1;		/* 0 = blasi_real */
										/* 1 = blasi_complex */

	unsigned int precision_type: 1;		/* 0 = blasi_single_precision */
										/* 1 = blasi_double_precision */

	unsigned int optimization : 3;		/* 0 = blasi_irregular  */
										/* 1 = blasi_regular   */
										/* 2 = blasi_block_irregular */
										/* 3 = blasi_block_regular */
										/* 4 = blasi_unassembled */
} SPBLASI_Matrix_prop_type;



/*
 *  These properties are independent of the implementation for 
 *	SPBLASI_Matrix.
 */

enum SPBLASI_property_values
{ 
	blasi_zero_base  			= 0, 
	blasi_one_base 				= 1,

	blasi_rowmajor				= 0,
	blasi_colmajor				= 1,

	blasi_no_repeated_indices	= 0,	/* this ensures that there are no	*/
	blasi_repeated_indices		= 1,	/* repeated indices, allowing for	*/
										/* possible optimziations.          */

	blasi_irregular				= 0,	/* In this implementation, these are */
	blasi_regular				= 1,	/* largely ignored.  Could provide   */
	blasi_block_regular			= 2,	/* opportunities for optimizations   */
	blasi_block_irregular		= 3,	/* in future versions.               */
	blasi_block_unassembled		= 4,

	blasi_general				= 0,	/* No special symmetry properties.   */
	blasi_lower_symmetric		= 1,	
	blasi_upper_symmetric		= 2,	
	blasi_lower_triangular		= 3,
	blasi_upper_triangular		= 4,
	blasi_lower_hermitian		= 5,
	blasi_upper_hermitian		= 6,

	blasi_real					= 0,
	blasi_complex				= 1,

	blasi_single_precision		= 0,
	blasi_double_precision		= 1
};


/**
* macros for SPBLASI_Matrix_property_tag.  The "get" version are pass
* by value, but the "set" versions are pass by reference, i.e.
*
* SPBLASI_is_zero_base(p);
* SPBLIAS_set_zero_base(&p);
*
*/

#define SPBLASI_is_zero_base(p)  ((p).index_base == blasi_zero_base)
#define SPBLASI_is_one_base(p)  ((p).index_base == blasi_one_base)

#define SPBLASI_is_rowmajor(p)  ((p).order == blasi_rowmajor)
#define SPBLASI_is_colmajor(p)  ((p).order == blasi_colmajor)

#define SPBLASI_has_repeated_indices(p)  \
	((p).repeated_indices == blasi_repeated_indices)


#define SPBLASI_is_double_precision(p) \
	((p).scalar_type == blasi_double_precision)
#define SPBLASI_is_complex(p)		 ((p).scalar_type == blasi_complex)
#define SPBLASI_is_real(p)	((p).scalar_type == blasi_real)

#define SPBLASI_is_general(p) ((p).structure == blasi_general)
#define SPBLASI_is_lower_symmetric(p)	((p).structure == blasi_lower_symmetric)
#define SPBLASI_is_upper_symmetric(p)	((p).structure == blasi_upper_symmetric)
#define SPBLASI_is_symmetric(p)	(SPBLASI_is_upper_symmetric(p) || \
							     SPBLASI_is_lower_symmetric(p) )
#define SPBLASI_is_lower_triangular(p)((p).structure == blasi_lower_triangular)
#define SPBLASI_is_upper_triangular(p) ((p).structure == blasi_upper_triangular)
#define SPBLASI_is_triangular(p)	(SPBLASI_is_upper_triangular(p) || \
							        SPBLASI_is_lower_triangular(p) )
#define SPBLASI_is_upper_hermitian(p)	((p).structure == blasi_upper_hermitian)
#define SPBLASI_is_lower_hermitian(p) ((p).structure == blasi_lower_hermitian)
#define SPBLASI_is_hermitian(p)	(SPBLASI_is_upper_hermitian(p) || \
							     SPBLASI_is_lower_hermitian(p) )

#define SPBLASI_set_general(p) ((p)->structure = blasi_general)
#define SPBLASI_set_lower_symmetric(p)	((p)->structure=blasi_lower_symmetric)
#define SPBLASI_set_upper_symmetric(p)	((p)->structure=blasi_upper_symmetric)
#define SPBLASI_set_lower_triangular(p) ((p)->structure=blasi_lower_triangular)
#define SPBLASI_set_upper_triangular(p) ((p)->structure=blasi_upper_triangular)
#define SPBLASI_set_lower_hermitian(p)	((p)->structure=blasi_lower_hermitian)
#define SPBLASI_set_upper_hermitian(p)	((p)->structure=blasi_upper_hermitian)


#define SPBLASI_set_zero_base(p)	((p)->index_base = blasi_zero_base)
#define SPBLASI_set_one_base(p)	((p)->index_base = blasi_one_base)
#define SPBLASI_set_real(p)		( (p)->scalar_type = blasi_real )
#define SPBLASI_set_complex(p)		( (p)->scalar_type = blasi_complex )
#define SPBLASI_set_double_precision(p) \
	((p)->precision_type = blasi_double_precision)
#define SPBLASI_set_single_precision(p) \
	((p)->precision_type = blasi_single_precision)

/* and so on ... */


/* view the current properties on one line, without newline */

void SPBLASI_Matrix_prop_print(const SPBLASI_Matrix_prop_type *P);

#endif


