/*****************************************************************************
 * 
 * File:    volume_update_kernels.cu
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: volume_update_kernels.cu 692 2011-09-13 07:18:19Z astivala $
 *
 * CUDA kernels to update the volume vector.
 *
 ****************************************************************************/

#include "volume_update_kernels.h"
#include "atomic_add_double_kernel.h"

/*
 * init_volume_kernel() - set volume vector to all zero
 *
 * Parmeters:
 *       link_volumes- vector of volumes on each link (dim num_edges)
 *       num_edges   - number of linkes 
 *
 */
__global__ void init_volume_kernel(double link_volumes[], long num_edges)
{
  // each thread does as many iterations as necessary to cover all edges
  // (usually we would want each thread to only do a single edgs but
  // this way any number of edges can be handled with any number of threads)
  const long tid = blockDim.x * blockIdx.x + threadIdx.x; // thread id
  for (long k =  tid; k < num_edges; k += gridDim.x * blockDim.x)
    link_volumes[k] = 0.0;
}



/*
 * link_volume_update_kernel() - increment volume on each edge accord to OD matrix
 *
 * Parameters:
 *   link_init_nodes - init node for each link k from net link data 
 *   num_start_nodes - number of start nodes (=zones)
 *   od_matrix  - demands array (OD matrix) in matrix form (square dim num_start_nodes)
 *                od_matrix[i * num_start_nodes  + j ] is deamdn from i to j
 *    link_volumes - (IN/OUT) - volume on each link
 *    predlink  -  predecessor link array,
 *
 * This kernel is run in paraell for each origin-destination pair and, given
 * the already computed cost/predecessor matrices, goes along the path
 * for the origin-destination pair and increments the volume for eadch
 * edge in the path by the demand for that O-D pair.
 *
 */

__global__ void link_volume_update_kernel(long link_init_nodes[],
                                          long num_start_nodes,
                                          double od_matrix[],
                                          double link_volumes[],
                                          long predlink[])
{
  long orig,dest,k,v;
  double routeflow;

  // each thread does as many iterations as necessary to cover all zones
  // (usually we would want each thread to only do a single zone but
  // this way any number of zones can be handled with any number of threads)
  for (orig = blockIdx.x+1; orig < num_start_nodes+1; orig += gridDim.x) 
  {
    for (dest = threadIdx.x+1; dest < num_start_nodes+1; dest += blockDim.x)
    {
      routeflow =  od_matrix[orig * (num_start_nodes+1) + dest];
      if (orig == dest || routeflow <= 0)
        continue;
      v = dest;
      while (v != orig)
      {
        k = predlink[v * num_start_nodes + orig];
        if (k < 0)
          break;  // should not happen; nonzero demand but no path
//        assert(k >= 0);
//        assert(links[k].term_node == v)
        atomicAdd(&link_volumes[k], routeflow);
        v = link_init_nodes[k];
      }
    }
  }
}
