#ifndef BLAS_MALLOC_H
#define BLAS_MALLOC_H 1

#include <stddef.h>
#include <malloc.h>


/* stddef is needed for size_t */
void  *blas_malloc(size_t size);
void *blas_realloc(void *p, size_t size);
void blas_free(void *p);


#endif
