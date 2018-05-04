#ifndef ATOMIC_ADD_DOUBLE_KERNEL_H
#define ATOMIC_ADD_DOUBLE_KERNEL_H
/*****************************************************************************
 * 
 * File:    atomic_add_double_kernel.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: atomic_add_double_kernel.h 682 2011-09-12 06:40:52Z astivala $
 *
 * CUDA atomic operation device functions.
 *
 * Requires device compute capabillity at least 1.2 (uses 64 bit atomic
 * functions on global memory).
 *
 ****************************************************************************/



/****************************************************************************
 * 
 * __device__ functions: callable on GPU only, inlined
 *
 ****************************************************************************/

// this function is from the NVIDIA Cuda Programming Guide (section B.11)
// we need it because there is no atomicAdd() for double precision floating point
// but it can be implemted with CAS on 64-bit word
__device__ double atomicAdd(double* address, double val)
{
  double old = *address, assumed;
  do {
    assumed = old;
    old =
      __longlong_as_double(
        atomicCAS((unsigned long long int*)address,
                  __double_as_longlong(assumed),
                  __double_as_longlong(val + assumed)));
  } while (assumed != old);
  return old;
}


#endif /* ATOMIC_ADD_DOUBLE_KERNEL_H */
