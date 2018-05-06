/*****************************************************************************
 * 
 * File:    okuyama_kernels.cu
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: okuyama_kernels.cu 224 2011-04-13 06:09:10Z astivala $
 *
 * CUDA kernels for multiple-source shortest path implementation using
 * CUDA based on:
 *
 * Okuyama, Ino, Hagihara 2008 "A Task Parallel Algorithm for Computing
 * the costs of All-Pairs Shortest Paths on the CUDA-compatible GPU"
 * Intl. Symp. Parallel Distributed Processing with Applications. 284-291
 *
 * Requires device compute capabillity at least 1.2 (uses 64 bit atomic
 * functions).
 * Developed on CUDA 2.3 (uses atomic functions).
 *
 * Each block is responsible for all start nodes (singe source
 * shortest path problems) for a single vertex. The Ma, Ca, Ua
 * (modification, cost, update) arrays are stored as Xa[i*num_nodes+j]
 * where j is source node and i is vertex (destination node). So
 * threads are blocked into all source nodes for one destination
 * vertex in a block and storage is arranged so all (destination)
 * vertices for a single source are contiguous. I.e. interleaving
 * for coalesced memory access.
 * 
 *
 ****************************************************************************/

#include <cutil_inline.h>      /* CUDA SDK */

#include "sssp.h"
#include "okuyama_kernels.h"
#include "atomic_devfunc.h"


/****************************************************************************
 * 
 * __constant__ memory: readonly memory from device
 *
 ****************************************************************************/

/* the Va array for packed adjancey list */
//__constant__ int c_Va[MAX_CONSTANT_NODES];

/****************************************************************************
 * 
 * __device__ data: global memory space on device
 *
 ****************************************************************************/

 // count updates from update kernel
__device__ unsigned int d_okuyama_update_count = 0;



/****************************************************************************
 * 
 * __global__ functions: GPU kernels, callable from host
 *
 ****************************************************************************/


/*
 * okuyama_init_mask_cost_update_arrays() - each thread for one vertex 
 *                           initialize mask, cost, and update arrays
 *
 * Parameters:
 *    Ma         - modification set: initialize to false
 *    Ca         - current cost:     intiialize to infinity 
 *    Ua         - updated cost and predecessor:     initialize to infinity,i
 *    start_nodes- source nodes:     Ma[s,s] = true, Ca[s,s] = 0
 *    num_nodes  - noumber of nodes
 *    num_start_nodes - number of source nodes
 *    Pa         - precessor node:   intialize Pa[i] = i
 *
 * Return value:
 *    None.
 */
__global__ void okuyama_init_mask_cost_update_arrays(bool Ma[], 
                                             float Ca[], 
                                                     cost_node_pair_t Ua[], 
                                                     int start_nodes[],
                                                     int num_nodes, 
                                                     int num_start_nodes,
                                                     int d_Pa[])
{
  // each thread does as many iterations as necessary to cover all nodes
  // (usually we would want each thread to only do a single node but
  // this way any number of nodes can be handled with any number of threads)
  for (int v = blockIdx.x; v < num_nodes; v += gridDim.x) 
  {
    for (int i = threadIdx.x; i < num_start_nodes; i += blockDim.x)
    {
      int s = start_nodes[i];
#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA)) 
      fprintf(stderr, "node %d source %d\n", v, s);
#endif
      if (v == s)
      {
        Ma[v * num_start_nodes + s] = true;
        Ca[v * num_start_nodes + s] = 0;
        Ua[v * num_start_nodes + s].cost = 0;
      }
      else
      {
        Ma[v * num_start_nodes + s] = false;
        Ca[v * num_start_nodes + s] = FLOATINF;
        Ua[v * num_start_nodes + s].cost = FLOATINF;
      }
      d_Pa[v * num_start_nodes + s] = v;
      d_Pa[v * num_start_nodes + s] = v;
    }
  }
}


/*
 * okuyama_scatter_kernel() - each thread for one vertex 
 *                           update costs of its adjcacent vertices 
 *
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list format
 *    Ma         - modification set
 *    Ca         - current cost
 *    Ua         - updated cost, predecessor
 *    start_nodes- source nodes
 *    num_nodes  - noumber of nodes 
 *    num_start_nodes - number of ousrce nodes
 *    Pa         - precessor node
 *
 * Return value:
 *    None.
 *
 * Sets global update count to 0 at end in preparation for update kernel next
 */
__global__ void  okuyama_scatter_kernel(int Va[],
                                       int Ea[],
                                       float Wa[],
                                       bool Ma[],
                                       float Ca[],
                                        volatile cost_node_pair_t Ua[], 
                                        int start_nodes[],
                                        int num_nodes,
                                        int  num_start_nodes,
                                        int Pa[])
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x; // thread id

  // TODO use shared memory - but maybe no use on Fermi architecture
  // anyway (already constrainted to use compute capability 2.0 for
  // atomic 64 bit etc) where caching means we get benefit of fast
  // memory without explicitly using shared memory?

  // each thread does as many iterations as necessary to cover all nodes
  // (usually we would want each thread to only do a single node but
  // this way any number of nodes can be handled with any number of threads)
  for (int v = blockIdx.x; v < num_nodes; v += gridDim.x) 
  {
    for (int k = threadIdx.x; k < num_start_nodes; k += blockDim.x)
    {
      int s = start_nodes[k];
      if (Ma[v*num_start_nodes+s])
      {
        Ma[v*num_start_nodes+s] = false;
        //for (int i = c_Va[v]; i < c_Va[v+1]; i++)  // all neighbours of v
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
          atomic_min_cost_node(&Ua[n*num_start_nodes+s],
                               Ca[v*num_start_nodes+s] + Wa[i], v);
        }
      }
    }
  }
  if (tid == 0) // only one thread needs to set global update count to 0
    d_okuyama_update_count = 0;
}

/*
 * okuyama_update_kernel() - each thread updates cost at its vertex
 *
 * Parameters:
 *    Ma         - modification set
 *    Ca         - current cost
 *    Ua         - updated cost, predecessor
 *    start_nodes- source nodes
 *    num_nodes  - noumber of nodes 
 *    num_start_nodes -= number of source nodes
 *    Pa         - precessor node
 *
 * Return value:
 *    None.
 *
 * Atomically updates the global updated count by atomic increment
 * if cost is modified. 
 */
__global__ void okuyama_update_kernel(
                                     bool Ma[], float Ca[], 
                                      cost_node_pair_t Ua[],
                                      int start_nodes[],
                                      int num_nodes,
                                      int num_start_nodes,
                                      int Pa[])
{
  // each thread does as many iterations as necessary to cover all nodes
  // (usually we would want each thread to only do a single node but
  // this way any number of nodes can be handled with any number of threads)
  for (int v = blockIdx.x; v < num_nodes; v += gridDim.x) 
  {
    for (int i = threadIdx.x; i < num_start_nodes; i += blockDim.x)
    {
      int s = start_nodes[i];
      if (Ca[v*num_start_nodes+s] > Ua[v*num_start_nodes+s].cost)
      {
        Ca[v*num_start_nodes+s] = Ua[v*num_start_nodes+s].cost;
        Ma[v*num_start_nodes+s] = true;
        Pa[v*num_start_nodes+s] = Ua[v*num_start_nodes+s].node;
        // don't understand atomicInc() on CUDA, using atomicAdd() instead
//XXX        atomicAdd(&d_okuyama_update_count, 1);
        d_okuyama_update_count = 1; // don't actually need count, just 0 or 1
      }
//XXX      Ua[v*num_start_nodes+s].cost = Ca[v*num_start_nodes+s];
//XXX      Ua[v*num_start_nodes+s].node = Pa[v*num_start_nodes+s];
    }
  }
}
