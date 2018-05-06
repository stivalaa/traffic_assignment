/*****************************************************************************
 * 
 * File:    harish_host.cu
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: harish_host.cu 224 2011-04-13 06:09:10Z astivala $
 *
 * CUDA host code for single-source shortest path implementation using
 * CUDA based on:
 *
 * Harish and Narayanan 2007 "Accelerating Large Graph Algorithms on the GPU
 * Using CUDA" HiPC 2007, LNCS 4873: 197-208
 * 
 *
 * Developed on CUDA 2.3.
 * Requires device compute capabillity at least 1.2 (uses 64 bit atomic
 * functions).
 *
 ****************************************************************************/

#include <assert.h>

#include <cutil_inline.h>      /* CUDA SDK */

#include "sssp.h"
#include "harish_kernels.h"
#include "harish_host.h"
#include "atomic_devfunc_types.h"

/*
 * harish_sssp() - single-source shortest path by Harish&Naryanan algorithm
 *
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    alloc_num_nodes - number of elements of Va allocated for (>=num_nodes)
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    start_node - source node
 *    distances (OUT) - array of shortest costs from source to each node
 *    predecessors (OUT) - arrya of predecessor nodes for each node
 *
 * Return value:
 *    None.
 *
 */
void harish_sssp(int Va[], int Ea[], float Wa[], 
                 int alloc_num_nodes,
                 int num_nodes, int num_edges,
                 int start_node, float distances[], int predecessors[])
{
  int *d_Va, *d_Ea;
  float *d_Wa;
  bool *d_Ma;
  float *d_Ca;
  volatile cost_node_pair_t *d_Ua;
  int *d_Pa;

  assert(sizeof(cost_node_pair_t) == 8); // must be exactly 64 bits

  dim3 dimBlock(num_nodes);    // threads per block
  dim3 dimGrid(1);             // blocks (per grid)
  if (num_nodes > BLOCK_SIZE) 
  {
    dimBlock = dim3(BLOCK_SIZE);
    dimGrid = dim3(alloc_num_nodes / BLOCK_SIZE);
  }

//  dimBlock = dim3(1); dimGrid = dim3(1); // XXX

//  fprintf(stderr, "Execution configuration: Grid = (%d,%d,%d) Block = (%d,%d,%d)\n", dimGrid.x,dimGrid.y,dimGrid.z, dimBlock.x,dimBlock.y,dimBlock.z);

  // allocate arrays for packed adjancey list format and 
  // copy graph in packed adjacney list format to device
  cutilSafeCall( cudaMalloc((void **)&d_Va, alloc_num_nodes*sizeof(int)) );
  cutilSafeCall( cudaMalloc((void **)&d_Ea, num_edges*sizeof(int)) );
  cutilSafeCall( cudaMalloc((void **)&d_Wa, num_edges*sizeof(float)) );

  cutilSafeCall( cudaMemcpy(d_Va, Va, alloc_num_nodes*sizeof(int),
                            cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy(d_Ea, Ea, num_edges*sizeof(int),
                            cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy(d_Wa, Wa, num_edges*sizeof(float),
                            cudaMemcpyHostToDevice) );

  // allocate arrays for modification set, cost, updated cost,
  // predcessor, updated predecessor
  cutilSafeCall( cudaMalloc((void **)&d_Ma, num_nodes*sizeof(bool)) );
  cutilSafeCall( cudaMalloc((void **)&d_Ca, num_nodes*sizeof(float)) );
  cutilSafeCall( cudaMalloc((void **)&d_Ua, num_nodes*sizeof(cost_node_pair_t)) );
  cutilSafeCall( cudaMalloc((void **)&d_Pa, num_nodes*sizeof(int)) );
  

  // initialize the  modification set, cost, updated cost arrays on device
  init_mask_cost_update_arrays<<<dimGrid, dimBlock>>>(d_Ma, d_Ca, 
                                                  (cost_node_pair_t *)d_Ua,
                                                      start_node,
                                                      d_Pa, 
                                                      num_nodes);
/*       cudaError_t cuda_errcode = cudaGetLastError(); */
/*       if (cuda_errcode != cudaSuccess) */
/*       { */
/*         fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cuda_errcode)); */
/*         exit(1); */

/*       } */

  CUT_CHECK_ERROR("Kernel execution failed (init_mask_cost_update_arrays)");
  cutilSafeCall( cudaThreadSynchronize() );

  // execute scatter kernel followed by update kernel while modified nodes
  unsigned int update_count = 1;
  unsigned int iter_count = 0;
  do
  {
    harish_scatter_kernel<<<dimGrid, dimBlock>>>(d_Va, d_Ea, d_Wa,
                                                 d_Ma, d_Ca, d_Ua,
                                                 d_Pa, num_nodes);
    CUT_CHECK_ERROR("Kernel execution failed (harish_scatter_kernel)");
    cutilSafeCall( cudaThreadSynchronize() );

    harish_update_kernel<<<dimGrid, dimBlock>>>(d_Va, d_Ea, d_Wa,
                                                d_Ma, d_Ca, 
                                                (cost_node_pair_t *)d_Ua,
                                                d_Pa, num_nodes);
    CUT_CHECK_ERROR("Kernel execution failed (harish_update_kernel)");
    cutilSafeCall( cudaThreadSynchronize() );

    cudaMemcpyFromSymbol(&update_count,"d_update_count",sizeof(unsigned int));

#ifdef DEBUG
    fprintf(stderr, "update_count = %d\n", update_count);
    printf("iter_count = %d\n", iter_count);
    cutilSafeCall( cudaMemcpy(distances, d_Ca, num_nodes*sizeof(float),
                            cudaMemcpyDeviceToHost) );
    for (int i = 0; i < num_nodes; i++)
      printf("%f ", distances[i]);
    printf("\n");
#endif /* DEBUG */

    iter_count++;
  }
  while (update_count > 0);

  // get the final costs and predcessors nodes back from the device
  cutilSafeCall( cudaMemcpy(distances, d_Ca, num_nodes*sizeof(float),
                            cudaMemcpyDeviceToHost) );
  cutilSafeCall( cudaMemcpy(predecessors, d_Pa, num_nodes*sizeof(int),
                            cudaMemcpyDeviceToHost) );

  // free device memory
  cutilSafeCall( cudaFree(d_Va) );
  cutilSafeCall( cudaFree(d_Ea) );
  cutilSafeCall( cudaFree(d_Wa) );
  cutilSafeCall( cudaFree(d_Ca) );
  cutilSafeCall( cudaFree((void *)d_Ua) );
  cutilSafeCall( cudaFree(d_Ma) );
  cutilSafeCall( cudaFree(d_Pa) );
}
