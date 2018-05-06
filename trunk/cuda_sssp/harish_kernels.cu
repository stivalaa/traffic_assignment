/*****************************************************************************
 * 
 * File:    harish_kernels.cu
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: harish_kernels.cu 224 2011-04-13 06:09:10Z astivala $
 *
 * CUDA kernels for single-source shortest path implementation using
 * CUDA based on:
 *
 * Harish and Narayanan 2007 "Accelerating Large Graph Algorithms on the GPU
 * Using CUDA" HiPC 2007, LNCS 4873: 197-208
 * 
 * Developed on CUDA 2.3.
 * Requires device compute capabillity at least 1.2 (uses 64 bit atomic
 * functions).
 *
 ****************************************************************************/

#include <cutil_inline.h>      /* CUDA SDK */

#include "sssp.h"
#include "harish_kernels.h"
#include "atomic_devfunc_types.h"
#include "atomic_devfunc.h"


/****************************************************************************
 * 
 * __device__ data: global memory space on device
 *
 ****************************************************************************/

__device__ unsigned int d_update_count = 0; // count updates from update kernel



/****************************************************************************
 * 
 * __global__ functions: GPU kernels, callable from host
 *
 ****************************************************************************/


/*
 * init_mask_cost_update_arrays() - each thread for one vertex 
 *                           initialize mask, cost, and update arrays
 *
 * Parameters:
 *    Ma         - modification set: initialize to false, except node s
 *    Ca         - current cost:     intiialize to infinity
 *    Ua         - updated cost and node:   initialize to infinity, i
 *    s          - start node:       Ma[s] = true, Ca[s] = 0
 *    Pa         - precessor node:   intialize Pa[i] = i
 *    num_nodes  - noumber of nodes (elements in Ma, Ca, Ua)
 *
 * Return value:
 *    None.
 */
__global__ void init_mask_cost_update_arrays(bool Ma[], 
                                             float Ca[], 
                                             cost_node_pair_t Ua[], int s,
                                             int Pa[],
                                             int num_nodes)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x; // thread id
  const int num_threads = blockDim.x * gridDim.x;  // total number of threads

  // each thread does as many iterations as necessary to cover all nodes
  // (usually we would want each thread to only do a single node but
  // this way any number of nodes can be handled with any number of threads)
  for (int v = tid; v < num_nodes; v += num_threads) 
  {
    if (v == s)
    {
      Ma[v] = true;
      Ca[v] = 0;
      Ua[v].cost = 0;
    }
    else
    {
      Ma[v] = false;
      Ca[v] = FLOATINF;
      Ua[v].cost = FLOATINF;
    }
    Pa[v] = v;
    Ua[v].node = v;
  }
}


/*
 * harish_scatter_kernel() - each thread for one vertex 
 *                           update costs of its adjcacent vertices 
 *
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list format
 *    Ma         - modification set
 *    Ca         - current cost
 *    Ua         - updated cost and node
 *    Pa         - precessor node
 *    num_nodes  - noumber of nodes (elements in Ma, Ca, Ua)
 *
 * Return value:
 *    None.
 *
 * Sets global update count to 0 at end in preparation for update kernel next
 */
__global__ void  harish_scatter_kernel(int Va[],
                                       int Ea[],
                                       float Wa[],
                                       bool Ma[],
                                       float Ca[],
                                       volatile cost_node_pair_t Ua[], 
                                       int Pa[],
                                       int num_nodes)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x; // thread id
  const int num_threads = blockDim.x * gridDim.x;  // total number of threads

  // each thread does as many iterations as necessary to cover all nodes
  // (usually we would want each thread to only do a single node but
  // this way any number of nodes can be handled with any number of threads)
  for (int v = tid; v < num_nodes; v += num_threads) 
  {
    if (Ma[v])
    {
      Ma[v] = false;
      for (int i = Va[v]; i < Va[v+1]; i++)  // all neighbours of v
      {
        int n = Ea[i];  // n is adjacent to v
#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA)) 
        fprintf(stderr, "node %d neighbour %d weight %f\n", v, n, Wa[i]);
#endif
        // As per Martin et al 2009 and Okuyama et al 2008 this must be atomic
        // however there is no floating point atomic min (atomic operations
        // on floating point are a problem, especially add etc. since
        // some commutative operators ar enot actually commutative on f.p
        // etc. ) so se use our own.
        // Also we need to update cost and predecessor node in a single
        // atomic instruction.
/*         if (Ca[v]+Wa[i] < Ua[n].cost) */
/*         { */
/*           Ua[n].cost = Ca[v] + Wa[i]; */
/*           Ua[n].node = v; */
/*         } */
        atomic_min_cost_node(&Ua[n], Ca[v] + Wa[i], v);
      }
    }
  }
  if (tid == 0) // only one thread needs to set global update count to 0
    d_update_count = 0;
}

/*
 * harish_update_kernel() - each thread updates cost at its vertex
 *
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list format
 *    Ma         - modification set
 *    Ca         - current cost
 *    Ua         - updated cost and node
 *    Pa         - precessor node
 *    num_nodes  - noumber of nodes (elements in Ma, Ca, Ua)
 *
 * Return value:
 *    None.
 *
 * Atomically updates the global updated count by atomic increment
 * if cost is modified. 
 */
__global__ void harish_update_kernel(int Va[], int Ea[], float Wa[],
                                     bool Ma[], float Ca[],
                                     cost_node_pair_t Ua[],
                                     int Pa[],
                                     int num_nodes)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x; // thread id
  const int num_threads = blockDim.x * gridDim.x;  // total number of threads
  
  // each thread does as many iterations as necessary to cover all nodes
  // (usually we would want each thread to only do a single node but
  // this way any number of nodes can be handled with any number of threads)
  for (int v = tid; v < num_nodes; v += num_threads) 
  {
    if (Ca[v] > Ua[v].cost)
    {
      Ca[v] = Ua[v].cost;
      Ma[v] = true;
      Pa[v] = Ua[v].node;
      // don't understand atomicInc() on CUDA, using atomicAdd() instead
//      atomicAdd(&d_update_count, 1);
      d_update_count = 1; // don't actually need count, just 0 or 1
    }
    Ua[v].cost = Ca[v];
    Ua[v].node = Pa[v];
  }
}
