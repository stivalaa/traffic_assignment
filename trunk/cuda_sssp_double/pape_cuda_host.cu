/*****************************************************************************
 * 
 * File:    pape_cuda_host.cu
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: pape_cuda_host.cu 781 2011-10-05 00:13:17Z astivala $
 *
 * CUDA host code for CUDA implemnetatnion of d'Esopo-Pape algorithm.
 *
 ****************************************************************************/

#include <assert.h>

#include <cutil_inline.h>      /* CUDA SDK */

#include "sssp.h"
#include "pape_kernels.h"
#include "pape_cuda_host.h"

extern bool using_fermi_architecture;

#define TIMER_DEBUG

//ceil(a / b)
extern "C" long iDivUp(int a, long b){
  return ((a % b) != 0) ? (a / b + 1) : (a / b);
}

/*
 * pape_cuda() - multiple-source shortest path by d'Esopo-Pape algorhtm.
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
 *  uses globals: using_fermi_architecture (readonly)
 *
 * Each CUDA thread does one source node.
 *
 */
void pape_cuda(long Va[], long Ea[], double Wa[], 
               long num_nodes, long num_edges,
               long start_nodes[], long num_start_nodes,
               double *distances, long *predecessors)
{
  long *d_Va = NULL, *d_Ea;
  double *d_Wa;
  long *d_start_nodes;
  long *d_Pa;
  double *d_Ca;
  unsigned int hTimer;
  double copytime, runtime;
  long first_thru_node =  0;
  long *d_queue_next;

  // allocate arrays for packed adjancey list format and 
  // copy graph in packed adjacney list format to device
  // also start nodes list
  cutilSafeCall( cudaMalloc((void **)&d_Va, (num_nodes+1)*sizeof(long)) );
  cutilSafeCall( cudaMalloc((void **)&d_Ea, num_edges*sizeof(long)) );
  cutilSafeCall( cudaMalloc((void **)&d_Wa, num_edges*sizeof(double)) );
  cutilSafeCall( cudaMalloc((void **)&d_start_nodes, num_start_nodes*sizeof(long)) );

  printf("%d nodes (%d KB) %d edges (%d KB)\n",
         num_nodes, 
         ( (num_nodes+1)*sizeof(long) ) / 1024,
         num_edges,
         ( num_edges*sizeof(long) + num_edges*sizeof(double) ) / 1024);

  cutilCheckError( cutCreateTimer(&hTimer) );
  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );

  cutilSafeCall( cudaMemcpy(d_Va, Va, (num_nodes+1)*sizeof(long),
                            cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy(d_Ea, Ea, num_edges*sizeof(long),
                            cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy(d_Wa, Wa, num_edges*sizeof(double),
                            cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy(d_start_nodes, start_nodes, num_start_nodes*sizeof(long),
                            cudaMemcpyHostToDevice) );
  
  cutilCheckError( cutStopTimer(hTimer) );
  copytime = cutGetTimerValue(hTimer);
  printf("time to copy %d nodes %d edges (total %d KB) to device: %f ms\n",
         num_nodes, num_edges,
         ( (num_nodes+1)*sizeof(long) + num_edges*sizeof(long) +
           num_edges*sizeof(double) ) / 1024,
         copytime);

  // allocate arrays for costs, predecessors, and queues
  cutilSafeCall( cudaMalloc((void **)&d_Ca, num_nodes*num_start_nodes*
                            sizeof(double)) );
  cutilSafeCall( cudaMalloc((void**)&d_Pa, num_nodes*num_start_nodes*
                            sizeof(long)) );
  cutilSafeCall( cudaMalloc((void **)&d_queue_next, 
                            num_nodes*num_start_nodes*sizeof(long)) );

  // initialize the  modification set, cost, updated cost arrays on device
  dim3 dimBlock(num_start_nodes);    // threads per block
  dim3 dimGrid(num_nodes);           // blocks (per grid)
  if (num_start_nodes > 512) // FIXME some rule for this
    dimBlock = dim3(512); // kernel will handle multiple start nodes per thread
   if (num_nodes > 65535) // FIXME some rule for this - 65535 is Fermi max */
     dimGrid = dim3(65535);
  fprintf(stdout, "Initalize Execution configuration: Grid = (%d,%d,%d) Block = (%d,%d,%d)\n", dimGrid.x,dimGrid.y,dimGrid.z, dimBlock.x,dimBlock.y,dimBlock.z);
#ifdef TIMER_DEBUG
    unsigned int tdhTimer;
    cutilCheckError( cutCreateTimer(&tdhTimer) );
    cutilCheckError( cutResetTimer(tdhTimer) );
    cutilCheckError( cutStartTimer(tdhTimer) );
#endif /* TIMER_DEBUG */

  pape_init_arrays<<<dimGrid, dimBlock>>>(num_nodes, num_start_nodes,
                                          d_start_nodes, d_Ca, d_Pa,
                                          d_queue_next);
  CUT_CHECK_ERROR("Kernel execution failed (okuyama_init_mask_cost_update_arrays)");
  cutilSafeCall( cudaThreadSynchronize() );

#ifdef TIMER_DEBUG
    cutilCheckError( cutStopTimer(tdhTimer) );
    double init_time = cutGetTimerValue(tdhTimer);
    fprintf(stderr, "pape_init_array time: %f ms\n", init_time);
#endif /* TIMER_DEBUG */



  dimBlock = dim3(num_start_nodes);    // threads per block
  dimGrid = dim3(1);           // blocks (per grid)
  if (num_start_nodes > 512) // FIXME some rule for this
  {
    dimBlock = dim3(512); // kernel will handle multiple start nodes per thread
    dimGrid = dim3(iDivUp(num_start_nodes, 512));
  }
  if (using_fermi_architecture)
  {
    // on Fermi (Tesla C2050) seem to get kernel failures or some kind of deadlock
    // with dimGrid < 16  (FIXME why?) - does not happn on GTX285, GTX280, Tesla C1060
    dimGrid = dim3(dimGrid.x < 16 ? 16 : dimGrid.x); 
  }
  fprintf(stdout, "Execution configuration: Grid = (%d,%d,%d) Block = (%d,%d,%d)\n", dimGrid.x,dimGrid.y,dimGrid.z, dimBlock.x,dimBlock.y,dimBlock.z);

#ifdef TIMER_DEBUG
    cutilCheckError( cutCreateTimer(&tdhTimer) );
    cutilCheckError( cutResetTimer(tdhTimer) );
    cutilCheckError( cutStartTimer(tdhTimer) );
#endif /* TIMER_DEBUG */
    
    pape_kernel<<<dimGrid, dimBlock>>>(d_Va, d_Ea, d_Wa,
                                       num_nodes, num_edges, num_start_nodes,
                                       d_start_nodes, first_thru_node,
                                       d_Ca, d_Pa, d_queue_next);

    CUT_CHECK_ERROR("Kernel execution failed (pape_kernel)");
    cutilSafeCall( cudaThreadSynchronize() );
#ifdef TIMER_DEBUG
    cutilCheckError( cutStopTimer(tdhTimer) );
    double pape_time = cutGetTimerValue(tdhTimer);
    fprintf(stderr, "pape_kernel time: %f ms\n", pape_time);
#endif /* TIMER_DEBUG */

  // get the final costs and predecessor nodes back from the device
  cutilSafeCall( cudaMemcpy(distances, d_Ca, 
                            num_nodes*num_start_nodes*sizeof(double),
                            cudaMemcpyDeviceToHost) );
  cutilSafeCall( cudaMemcpy(predecessors, d_Pa,
                            num_nodes*num_start_nodes*sizeof(long),
                            cudaMemcpyDeviceToHost) );

  // free device memory
  cutilSafeCall( cudaFree(d_Ca) );
  cutilSafeCall( cudaFree(d_Pa) );
  cutilSafeCall( cudaFree(d_start_nodes) );
  cutilSafeCall( cudaFree(d_Va) );
  cutilSafeCall( cudaFree(d_Ea) );
  cutilSafeCall( cudaFree(d_Wa) );
  cutilSafeCall( cudaFree(d_queue_next) );

}
