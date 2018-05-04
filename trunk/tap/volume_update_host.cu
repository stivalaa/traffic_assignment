/*****************************************************************************
 * 
 * File:    volume_update_host.cu
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: volume_update_host.cu 692 2011-09-13 07:18:19Z astivala $
 *
 * CUDA host code for volume vector update
 *
 ****************************************************************************/

#include <assert.h>

#include <cutil_inline.h>      /* CUDA SDK */

#include "volume_update_kernels.h"
#include "volume_update_host.h"
#include "utils.h"

#define TIMER_DEBUG
#define VERBOSE

/*****************************************************************************
 * 
 * static data
 *
 ****************************************************************************/

// links and OD matrix device data set in link_volumne_data_setup()
static long *d_link_init_nodes = NULL;
static double *d_od_matrix = NULL;

/*****************************************************************************
 * 
 * functions
 *
 ****************************************************************************/

/*
 * link_volume_data_setup() - set the link and O-D demand matrix data on device
 *
 * Parameters:
 *   links - net link data 
 *    num_links - number of links
 *   num_start_nodes - number of start nodes (=zones)
 *    demands  - demands array parsed from trips file
 *                for each origin demands[origin] is an array of demands structs
 *                terminated by one with 0  dest (origins and dests start at 1)
 *
 *  Return value;
 *    None.
 */

void link_volume_data_setup(link_data_t *links, long num_links,
                            long num_start_nodes,
                            demand_data_t *demands[])
{
  // convert dmeands data to 2d matrix (no pointers)
  double *od_matrix = (double *)calloc((num_start_nodes+1)*(num_start_nodes+1), sizeof(double));
  if (!od_matrix)
  {
    fprintf(stderr, "calloc od_matrix failed\n");
    exit(1);
  }
  long   dest;
  double routeflow;
  for (long orig = 1; orig < num_start_nodes+1; orig++)
  {
    for (long i = 0; (dest = demands[orig][i].dest) != 0; i++)
    {
      assert(dest > 0);
      assert(dest <= num_start_nodes);
      routeflow = demands[orig][i].demand;
      od_matrix[orig * (num_start_nodes+1) + dest] = routeflow;
    }
  }

  // build vector of link init_nodes, we only need that on device not all link data
  long *link_init_nodes = (long *)malloc(num_links*sizeof(long));
  if (!link_init_nodes) 
  {
    fprintf(stderr, "malloc link_init_nodes failed\n");
    exit(1);
  }
  for (long k = 0; k < num_links; k++)
    link_init_nodes[k] = links[k].init_node;

  // allocate device memory and send matrix and vector to device

  cutilSafeCall( cudaMalloc((void **)&d_link_init_nodes, num_links*sizeof(long)) );
  cutilSafeCall( cudaMalloc((void **)&d_od_matrix, (num_start_nodes+1)*(num_start_nodes+1)*sizeof(double)) );

#ifdef TIMER_DEBUG
  fprintf(stderr, "%d x %d OD matrix (%d KB)\n", num_start_nodes+1, 
          num_start_nodes+1,
          (num_start_nodes+1)*(num_start_nodes+1)*sizeof(double)/1024);
  fprintf(stderr, "link init_nodes for %d links (%d KB)\n", num_links,
          num_links*sizeof(long)/1024);
  unsigned int hTimer;
  double copytime;
  cutilCheckError( cutCreateTimer(&hTimer) );
  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );
#endif

  cutilSafeCall( cudaMemcpy(d_od_matrix, od_matrix, 
                            (num_start_nodes +1)* (num_start_nodes+1) * sizeof(double),
                            cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy(d_link_init_nodes, link_init_nodes, num_links*sizeof(long),
                            cudaMemcpyHostToDevice) );

#ifdef TIMER_DEBUG
  cutilCheckError( cutStopTimer(hTimer) );
  copytime = cutGetTimerValue(hTimer);
  fprintf(stderr, "time to copy OD matrix and link data to device: %f ms\n", 
          copytime);
#endif                            
  
  free(od_matrix);
  free(link_init_nodes);
}


/*
 * link_volume_data_cleanup() - sfree link and O-D demand matrix dataon device
 *
 * Parameters:
 *    None.
 *  Return value;
 *    None.
 */
void link_volume_data_cleanup(void)
{
  cutilSafeCall( cudaFree(d_link_init_nodes) );
  cutilSafeCall( cudaFree(d_od_matrix) );
}


/*
 * link_volume_update - update the link volume vector CUDA host code
 *
 * Using the results from the shortest path computatinos, update the volume
 * on each path. This version is CUDA host code for running kernel on GPU.
 *
 * The links and demand data must already be setup on GPU (by init_linkvolume_data())
 *
 *   num_start_nodes - number of start nodes (=zones)
 *   num_nodes = number of nodes
 *    num_edges - number of links (=edges)
 *    link_volumes - (OUT) - volume on each link
 *    predlink  -  predecessor link array,
 *
 *  Return value;
 *    None.
 *
 */

void link_volume_update(long num_start_nodes, long num_edges, long num_nodes,
                        double link_volumes[],
                        long predlink[])
{

 // TODO FIXME should use the predlink already on device here instead of copying back again

  // allocate arrays on device
  long *d_predlink = NULL;
  double *d_link_volumes = NULL;
  cutilSafeCall( cudaMalloc((void **)&d_predlink,
                            (num_nodes+1)*num_start_nodes*sizeof(long)) );
  cutilSafeCall( cudaMalloc((void **)&d_link_volumes,
                            num_edges * sizeof(double)) );

  // copy arrays to device
  cutilSafeCall( cudaMemcpy(d_predlink, predlink, 
                            (num_nodes+1)*num_start_nodes*sizeof(long),
                            cudaMemcpyHostToDevice) );
  
  dim3 dimBlock, dimGrid;

  if (num_edges < 512)
  {
    dimBlock = dim3(num_start_nodes);
    dimGrid = dim3(1);
  }
  else
  {
    dimBlock = dim3(512); // threads per block
    dimGrid = dim3(iDivUp(num_edges, 512)); // blocks (per grid)
  }

#ifdef VERBOSE
  fprintf(stdout, "init_volume_kernel Execution configuration: Grid = (%d,%d,%d) Block = (%d,%d,%d)\n", dimGrid.x,dimGrid.y,dimGrid.z, dimBlock.x,dimBlock.y,dimBlock.z);
#endif
#ifdef TIMER_DEBUG
  unsigned int hTimer;
  cutilCheckError( cutCreateTimer(&hTimer) );
  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );
#endif

 init_volume_kernel<<<dimGrid, dimBlock>>>(d_link_volumes, num_edges);
 CUT_CHECK_ERROR("kernel executtino failed (init_volume_kernel)");
 cutilSafeCall( cudaThreadSynchronize() );

#ifdef TIMER_DEBUG
  cutilCheckError( cutStopTimer(hTimer) );
  double xtime = cutGetTimerValue(hTimer);
  fprintf(stdout,  "init link volume time %f ms\n", xtime);
#endif


 dimBlock = dim3(num_start_nodes); // threads per block
 dimGrid = dim3(num_start_nodes) ; // blocks (per grid)
 if (num_start_nodes > 128) // FIXME some rule for this
   dimBlock = dim3(128);


#ifdef VERBOSE
  fprintf(stdout, "link_volume_update_kernel Execution configuration: Grid = (%d,%d,%d) Block = (%d,%d,%d)\n", dimGrid.x,dimGrid.y,dimGrid.z, dimBlock.x,dimBlock.y,dimBlock.z);
#endif
#ifdef TIMER_DEBUG
  cutilCheckError( cutCreateTimer(&hTimer) );
  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );
#endif
 link_volume_update_kernel<<<dimGrid, dimBlock>>>(d_link_init_nodes, num_start_nodes, d_od_matrix,
                                                  d_link_volumes, d_predlink);
 CUT_CHECK_ERROR("kernel executtino failed (link_volume_update_kernel)");
 cutilSafeCall( cudaThreadSynchronize() );
 
#ifdef TIMER_DEBUG
  cutilCheckError( cutStopTimer(hTimer) );
  xtime = cutGetTimerValue(hTimer);
  fprintf(stdout,  "link volume update time %f ms\n", xtime);
#endif


 // get results back from device
 cutilSafeCall( cudaMemcpy( link_volumes, d_link_volumes, 
                            num_edges * sizeof(double),
                            cudaMemcpyDeviceToHost) );
                           
  // free device memory allocaed here
  cutilSafeCall( cudaFree(d_predlink) );
  cutilSafeCall( cudaFree(d_link_volumes) );
}
