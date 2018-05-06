#include <stdio.h>
#include "blas_malloc.h"
#include "spblasi_matrix.h"
#include "table.h"


/* -------------- BEGIN OF LOCAL VARIABLE DECLARATIONS --------------------*/

/*
* This defines the global table of sparse matrices from which we
* active during program execution.  The SPBLASI_table_size
* is also used as an init flag. If it is 0, then no 
* Sparse BLAS matrices are currently allocated.  That is,
* the Sparse BLAS framework is not consuming memory.
*
* Note that in a threaded environment, making changes to this
* variable needs to be in a critical section.  These
* are documented appropriately.
*
*/


static Table *SPBLASI_Table = NULL;

/* 
* This avoids the user having to call some init() function to
* build the global tables for the Sparse BLAS matrices.
* If no Level 2 or Level 3 primitives are needed, the matrix
* (handle) framework is not contructed.
*
* The number of handles currently active (i.e. created via new,
* but not yet deleted.)   This is also used as a flag, when
* it becomes greater than 0, memory for the handle-management
* framework is allocated.  When the number dwindles down to zero,
* all associated memory (for global table of pointers) is released.
*
*/





/* -------------- END OF LOCAL VARIABLE DECLARATIONS --------------------*/




int SPBLASI_table_size(void)
{
	return Table_size(SPBLASI_Table);
}

int SPBLASI_table_num_active_matrices(void)
{
	return Table_num_entries(SPBLASI_Table);
}

/*
* Returns the internal sparse matrix from global table, or NULL if
* not found.
*
*/
SPBLASI_Matrix* SPBLASI_table_get(int i)
{
	return (SPBLASI_Matrix *) Table_get(SPBLASI_Table, i);
}



int SPBLASI_table_insert(SPBLASI_Matrix *A)
{
	if (SPBLASI_Table == NULL)
		SPBLASI_Table = Table_new(1);

	return Table_insert(SPBLASI_Table, A);
}     

void SPBLASI_table_destroy(void)
{
	Table_delete(SPBLASI_Table);
	SPBLASI_Table = NULL;
}


/**
* Remove matrix from global list, but does NOT destroy it.
* If not found, or if index "i" is out of range, function
* returns NULL, otherwise it returns a pointer to the sparse
* matrix taken from that spot.
*
*/
SPBLASI_Matrix *SPBLASI_table_remove(int i)
{
	SPBLASI_Matrix *t = Table_remove(SPBLASI_Table, i);
	if (SPBLASI_table_num_active_matrices() < 1)
		SPBLASI_table_destroy();

	return t;
}


void SPBLASI_table_dump(void)
{
	int i=0;
	int SPBLASI_Table_size = Table_size(SPBLASI_Table);
	printf("Sparse BLAS Table Dump:  %d item(s)\n", SPBLASI_Table_size);


	for (i=0; i< SPBLASI_Table_size; i++)
	{
		SPBLASI_Matrix *S = SPBLASI_table_get(i); 	
		printf("Slot [%i]: ", i);
		if (S)
			SPBLASI_Matrix_print_info(S);
		else
			printf(" empty.");
		printf("\n");
	}
	printf("\n");

		
}


#if 0

/* ------------------------ DIAGNOSTIC UTILITIES ----------------------- */

#include <stdio.h>
void SPBLASID_Matrix_dump(SPBLASI_Matrix S, int debug_level)
{
	DCSR A = NULL;
	SPBLASI_Matrix_prop_struct prop;

	if (S == NULL)
	{
		fprintf(stderr, "NULL SPBLASI_Matrix.\n");
		return;
	}
	A = S->mat_.dcsr_;
	prop = S->prop_;

	if (debug_level < 2) return;

	if (A == NULL)
		fprintf(stderr, "DCSR is NULL.\n");
	else
		fprintf(stderr, "M = %d, N = %d, nz = %d\n", DCSR_M(A), DCSR_N(A),
							DCSR_nz(A));
	fprintf(stderr, "\t mtype = %d\n", S->mtype_);
	fprintf(stderr, "\t state = %d\n", S->state_);
	fprintf(stderr, "\t prop: [%d:%d:%d:%d:%d]\n",
			prop.index_base, prop.order, prop.repeated_indices,
			prop.optimization, prop.scalar_type);

	if (debug_level < 3) return;

	DCSR_print(A, "%4.2g");
	fprintf(stderr, "\n");

}
	


void SPBLASI_delete_Matrix(int h)
{
	SPBLASI_Matrix S = SPBLASI_get_Matrix(h);

	if (S==NULL) 
		return;

	if (S->state_ == unused)
		return;

	DCSR_free(S->mat_.dcsr_);
	S->state_ = unused;
	SPBLASI_Table_num_active_matrices--;

	if (SPBLASI_Table_num_active_matrices < 1)
	{
		blas_free(SPBLASI_Table);
		SPBLASI_Table_size = 0;
	}
}



#endif

