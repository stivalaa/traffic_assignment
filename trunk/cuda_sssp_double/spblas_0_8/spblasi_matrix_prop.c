#include <stdio.h>
#include "spblasi_matrix_prop.h"


void SPBLASI_Matrix_prop_print(const SPBLASI_Matrix_prop_type *P)
{
	
	switch (P->structure)
	{
		case blasi_general:			printf("[general] ");
									break;

		case blasi_lower_symmetric:	printf("[lower-symmetric] ");
									break;

		case blasi_upper_symmetric:	printf("[upper-symmetric] ");
									break;

		case blasi_lower_hermitian:	printf("[lower-hermitian] ");
									break;

		case blasi_upper_hermitian:	printf("[upper-hermitian] ");
									break;

		case blasi_lower_triangular: printf("[lower-triangular] ");
									 break;

		case blasi_upper_triangular: printf("[upper-triangular] ");
									 break;

		default:					 printf("[error: invalid structure] ");
									 return;
	}

	switch (P->scalar_type)
	{
		case blasi_real:			printf("[real] ");
									break;

		case blasi_complex:			printf("[complex] ");
									break;

		default:					printf("[error: invalid scalar type] ");
									return;
	}

	switch (P->index_base)
	{
		case 0:						printf("[0-based] ");
									break;

		case 1:						printf("[1-based] ");
									break;
	
		default:					printf("[error: invalid index base] ");
									return;
	}


}
