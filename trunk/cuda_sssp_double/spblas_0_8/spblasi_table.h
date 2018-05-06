
#ifndef SPBLASI_TABLE_H
#define SPBLASI_TABLE_H 1

#include "spblasi_matrix.h"

int                  SPBLASI_table_insert(SPBLASI_Matrix *A);
SPBLASI_Matrix   	*SPBLASI_table_get(int i);
SPBLASI_Matrix      *SPBLASI_table_remove(int i);
void				 SPBLASI_table_destroy(void);
int                  SPBLASI_table_size(void);
int					 SPBLASI_table_num_active_matrices(void);
void				 SPBLASI_table_dump(void);

#endif
