/*****************************************************************************
 * 
 * File:    okuyama_host.cu
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: okuyama_host.cu 224 2011-04-13 06:09:10Z astivala $
 *
 * CUDA host code for multiple-source shortest path implementation using
 * CUDA based on:
 *
 * Okuyama, Ino, Hagihara 2008 "A Task Parallel Algorithm for Computing
 * the costs of All-Pairs Shortest Paths on the CUDA-compatible GPU"
 * Intl. Symp. Parallel Distributed Processing with Applications. 284-291
 *
 * Developed on CUDA 2.3.
 * Requires device compute capabillity at least 1.2 (uses 64 bit atomic
 * functions).
 *
 ****************************************************************************/

#include <assert.h>

#include <cutil_inline.h>      /* CUDA SDK */

#include "sssp.h"
#include "okuyama_kernels.h"
#include "okuyama_host.h"

#undef TIMER_DEBUG

/*
 * okuyama_sssp() - multiple-source shortest path by Okuyama et al algorithm
 *
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    start_nodes - array of source nodes
 *    num_start_nodes - number of source nodes (elements in start_nodes)
 *    distances (OUT) - 2d array of shortest costs from sources to each node
 *                      distances[i*num_nodes+j] is cost from sourc j to node i
 *    predecessors (OUT) - 2d array (as above)
 *                        of predecessor nodes for each node
 *
 * Return value:
 *    None.
 *
 * Each block is responsible for all start nodes (singe source shortest path
 * problems) for a single vertex.
 *
 */
void okuyama_sssp(int Va[], int Ea[], float Wa[], 
                  int num_nodes, int num_edges,
                  int start_nodes[], int num_start_nodes,
                  float *distances, int *predecessors)
{
  int *d_Va = NULL, *d_Ea;
  float *d_Wa;
  bool *d_Ma;
  float *d_Ca;
  volatile cost_node_pair_t *d_Ua;
  int *d_start_nodes;
  int *d_Pa;
  unsigned int hTimer;
  double copytime, runtime;

 // assert(num_nodes + 1 < MAX_CONSTANT_NODES); // TODO handle this if too many
  // using constant memory for this makes little improvement anyway

  dim3 dimBlock(num_start_nodes);    // threads per block
  dim3 dimGrid(num_nodes);           // blocks (per grid)
  if (num_start_nodes > 512) // FIXME some rule for this
    dimBlock = dim3(512); // kernel will handle multiple start nodes per thread

   if (num_nodes > 65535) // FIXME some rule for this - 65535 is Fermi max */
     dimGrid = dim3(65535);

  fprintf(stdout, "Execution configuration: Grid = (%d,%d,%d) Block = (%d,%d,%d)\n", dimGrid.x,dimGrid.y,dimGrid.z, dimBlock.x,dimBlock.y,dimBlock.z);

  // allocate arrays for packed adjancey list format and 
  // copy graph in packed adjacney list format to device
  // also start nodes list
  cutilSafeCall( cudaMalloc((void **)&d_Va, (num_nodes+1)*sizeof(int)) );
  cutilSafeCall( cudaMalloc((void **)&d_Ea, num_edges*sizeof(int)) );
  cutilSafeCall( cudaMalloc((void **)&d_Wa, num_edges*sizeof(float)) );
  cutilSafeCall( cudaMalloc((void **)&d_start_nodes, num_start_nodes*sizeof(int)) );

  printf("%d nodes (%d KB) %d edges (%d KB)\n",
         num_nodes, 
         ( (num_nodes+1)*sizeof(int) ) / 1024,
         num_edges,
         ( num_edges*sizeof(int) + num_edges*sizeof(float) ) / 1024);

  cutilCheckError( cutCreateTimer(&hTimer) );
  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );

  cutilSafeCall( cudaMemcpy(d_Va, Va, (num_nodes+1)*sizeof(int),
                            cudaMemcpyHostToDevice) );
  // copy Va to constant memory
//  cutilSafeCall( cudaMemcpyToSymbol("c_Va", Va, (num_nodes+1)*sizeof(int)) ); 
  // and others to global memory
  cutilSafeCall( cudaMemcpy(d_Ea, Ea, num_edges*sizeof(int),
                            cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy(d_Wa, Wa, num_edges*sizeof(float),
                            cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy(d_start_nodes, start_nodes, num_start_nodes*sizeof(int),
                            cudaMemcpyHostToDevice) );
  
  cutilCheckError( cutStopTimer(hTimer) );
  copytime = cutGetTimerValue(hTimer);
  printf("time to copy %d nodes %d edges (total %d KB) to device: %f ms\n",
         num_nodes, num_edges,
         ( (num_nodes+1)*sizeof(int) + num_edges*sizeof(int) +
           num_edges*sizeof(float) ) / 1024,
         copytime);

  // allocate arrays for modification set, cost, updated cost
  //  and predecessor arrays
  cutilSafeCall( cudaMalloc((void **)&d_Ma, 
                            num_nodes*num_start_nodes*sizeof(bool)) );
  cutilSafeCall( cudaMalloc((void **)&d_Ca, num_nodes*num_start_nodes*
                            sizeof(float)) );
  cutilSafeCall( cudaMalloc((void **)&d_Ua, num_nodes*num_start_nodes*
                            sizeof(cost_node_pair_t)) );
  cutilSafeCall( cudaMalloc((void**)&d_Pa, num_nodes*num_start_nodes*
                            sizeof(int)) );

  // initialize the  modification set, cost, updated cost arrays on device
  okuyama_init_mask_cost_update_arrays<<<dimGrid, dimBlock>>>(d_Ma, d_Ca, 
                                                 (cost_node_pair_t *)d_Ua,
                                                              d_start_nodes, 
                                                              num_nodes,
                                                              num_start_nodes,
                                                              d_Pa);
  CUT_CHECK_ERROR("Kernel execution failed (okuyama_init_mask_cost_update_arrays)");
  cutilSafeCall( cudaThreadSynchronize() );

  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );

  // execute scatter kernel followed by update kernel while modified nodes
  unsigned int update_count = 1;
  unsigned int iter_count = 0;
  do
  {
#ifdef TIMER_DEBUG
    unsigned int tdhTimer;
    cutilCheckError( cutCreateTimer(&tdhTimer) );
    cutilCheckError( cutResetTimer(tdhTimer) );
    cutilCheckError( cutStartTimer(tdhTimer) );
#endif /* TIMER_DEBUG */
    okuyama_scatter_kernel<<<dimGrid, dimBlock>>>(d_Va,//NULL,//[use c_Va now] d_Va,
                                                  d_Ea, d_Wa,
                                                  d_Ma, d_Ca, d_Ua,
                                                  d_start_nodes,
                                                  num_nodes,
                                                  num_start_nodes,
                                                  d_Pa);
    CUT_CHECK_ERROR("Kernel execution failed (okuyama_scatter_kernel)");
    cutilSafeCall( cudaThreadSynchronize() );
#ifdef TIMER_DEBUG
    cutilCheckError( cutStopTimer(tdhTimer) );
    double scatter_time = cutGetTimerValue(tdhTimer);
    fprintf(stderr, "okuyama_scatter_kernel time: %f ms\n", scatter_time);
#endif /* TIMER_DEBUG */


#ifdef TIMER_DEBUG
    cutilCheckError( cutResetTimer(tdhTimer) );
    cutilCheckError( cutStartTimer(tdhTimer) );
#endif /* TIMER_DEBUG */
    okuyama_update_kernel<<<dimGrid, dimBlock>>>(
                                                 d_Ma, d_Ca, 
                                                 (cost_node_pair_t *)d_Ua,
                                                 d_start_nodes,
                                                 num_nodes,
                                                 num_start_nodes,
                                                 d_Pa);
    CUT_CHECK_ERROR("Kernel execution failed (okuyama_update_kernel)");
    cutilSafeCall( cudaThreadSynchronize() );
#ifdef TIMER_DEBUG
    cutilCheckError( cutStopTimer(tdhTimer) );
    double update_time = cutGetTimerValue(tdhTimer);
    fprintf(stderr, "okuyama_update_kernel time: %f ms\n", update_time);
#endif /* TIMER_DEBUG */

    cudaMemcpyFromSymbol(&update_count,"d_okuyama_update_count",
                         sizeof(unsigned int));

#ifdef DEBUG
    fprintf(stderr, "okuyama_update_count = %d\n", update_count);
    printf("iter_count = %d\n", iter_count);
#endif /* DEBUG */

    iter_count++;
  }
  while (update_count > 0);

  cutilCheckError( cutStopTimer(hTimer) );
  runtime = cutGetTimerValue(hTimer);
  printf("time to run %d iterations on device: %f ms\n",
         iter_count, runtime);

  // get the final costs and predecessor nodes back from the device
  cutilSafeCall( cudaMemcpy(distances, d_Ca, 
                            num_nodes*num_start_nodes*sizeof(float),
                            cudaMemcpyDeviceToHost) );
  cutilSafeCall( cudaMemcpy(predecessors, d_Pa,
                            num_nodes*num_start_nodes*sizeof(int),
                            cudaMemcpyDeviceToHost) );

  // free device memory
  cutilSafeCall( cudaFree(d_Pa) );
  cutilSafeCall( cudaFree(d_start_nodes) );
  cutilSafeCall( cudaFree(d_Va) );
  cutilSafeCall( cudaFree(d_Ea) );
  cutilSafeCall( cudaFree(d_Wa) );
  cutilSafeCall( cudaFree(d_Ca) );
  cutilSafeCall( cudaFree((void *)d_Ua) );
  cutilSafeCall( cudaFree(d_Ma) );


}
