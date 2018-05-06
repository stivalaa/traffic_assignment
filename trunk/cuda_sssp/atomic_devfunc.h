#ifndef ATOMIC_DEVFUNC_H
#define ATOMIC_DEVFUNC_H
/*****************************************************************************
 * 
 * File:    atomic_devfunc.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: okuyama_kernels.cu 90 2011-02-20 23:55:38Z astivala $
 *
 * CUDA atomic operation device functions.
 *
 * Requires device compute capabillity at least 1.2 (uses 64 bit atomic
 * functions).
 *
 * 
 * FIXME require -G0 or -G1 or these atomic functions don't work
 * (cause deadlock). Something to do with the way SIMT warps work?
 * If compiling with optimization (e.g. -O3 without -G ) it deadlocks
 * but -G0 or -G1 disables -O and it works in debug mode.
 *
 ****************************************************************************/

#include "atomic_devfunc_types.h"



/****************************************************************************
 * 
 * __device__ functions: callable on GPU only, inlined
 *
 ****************************************************************************/


/* 
 * atomic_min_float() - atomic min on floating point
 *
 * Parameters:
 *     address - address of float to do min on
 *     val     - value to compare to that at address
 *
 * Return value:
 *     Old value
 *
 *  Read float at address (global or shared memory) as old
 *  and store min(old,val) at address "atomically".
 *
 *  We need atomic min on floating point however there is no floating
 *  point atomic min (atomic operations on floating point are a
 *  problem, especially add etc. since some commutative operators ar
 *  enot actually commutative on f.p etc. ) so se use our own. We implement
 *  it using CAS logic on the float just treating it as 32 bit
 *  word.
 *
 */
__device__ float atomic_min_float(float *address, float val)
{
  float old;
  do
  {
    old = *address;
  }  
  while (__int_as_float(atomicCAS((unsigned int *)address, 
                                  __float_as_int(old),
                                  __float_as_int(fminf(old, val)))) != old);
  return old;
}



/* 
 * atomic_min_cost_node() - atomic min (cost, node) (float, int) pair.
 *
 * Parameters:
 *     address - address of (cost,node) pair to do min on
 *     val     - cost value to do min of cost on the pair
 *     newnode - new node value to set if val is the min
 *
 * Return value:
 *     Old value
 *
 *  Read float from (cost,node) pair at address (global or shared memory)
 *  as old and store (min(old,val), newnode) at address "atomically".
 *
 * Comparing equality of floating point numbers here by just comparing
 * their bit patters (as unsigned ints) could go very wrong with +0,-0
 * or denormalized numbers etc. though...?
 */
__device__ float atomic_min_cost_node(volatile cost_node_pair_t *address, 
                                      float val, int newnode)
{
  union __align__(8) costnode_longlong_union {
    cost_node_pair_t       cost_node_pair;  // (float,int)
    unsigned long long int ull; 
  };
  union costnode_longlong_union old;
  union costnode_longlong_union assumed;
  union costnode_longlong_union newvalue;
    
  old.cost_node_pair = *(cost_node_pair_t*)address;
  do   
  {
    assumed = old;
    newvalue.cost_node_pair.cost = fminf(old.cost_node_pair.cost, val);
    newvalue.cost_node_pair.node = val < old.cost_node_pair.cost 
                                                     ? newnode 
                                                     : old.cost_node_pair.node;
  }
  while ((old.ull = atomicCAS((unsigned long long int *)address,
                              assumed.ull, newvalue.ull))
         != assumed.ull);
  return old.cost_node_pair.cost;
}

#endif /* ATOMIC_DEVFUNC_H */
