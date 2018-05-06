/**
*
* Sparse BLAS error handler: this header file is implementation dependent
* and is NOT part of the standard specification.
*
* The default behavior is to return a non-zero value for BLAS requests
* that did not complete successfully.  As a debugging aid, this
* implementation allows one to redefine this behavior to use the
* ANSI C ASSERT macro to diagnose errors, by using the BLAS_ERRORS_ARE_FATAL
* flag at compile time.
*
* BLAS_ASSERT_RETURN(x, ret_val)  tests for condition "x" and it is FALSE,
* exits the function prematurely with a return value of "ret_val".
* However, if the flag BLAS_ERRORS_ARE_FATAL, then it will stop execution
* and print out a diagnostic message.
* 
* For example, 
*
* BLAS_ASSERT_RETURN(A != NULL, NULL);
*
* if the BLAS_ERRORS_ARE_FATAL flag is on, then the above is equivalent
* to
*
* assert(A != NULL);
*
* otherwise, it is equivalent to 
*
* if (A == NULL) return NULL;
*
*
*
*
*/


#ifndef SPBLASI_ERROR_H_
#define SPBLASI_ERROR_H_

#ifdef BLAS_ERRORS_ARE_FATAL
#include <assert.h>
#define BLAS_ASSERT_RETURN(x, ret_val) assert(x)
#else

#define BLAS_ASSERT_RETURN(x, ret_val) {if (!(x)) return ret_val;}
#endif



#endif
