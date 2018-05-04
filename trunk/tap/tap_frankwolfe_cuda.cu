/*****************************************************************************
 * 
 * File:    tap_frankwolfe_cuda.cu
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_frankwolfe_cuda.cu 823 2011-10-31 04:35:15Z astivala $
 *
 * Traffic assignment by Frank-Wolfe algorithm.
 *
 * There are many descriptions of this algorithm, a simple one for
 * example is in Sheffi, Y. 1985 "Urban Transportation Networks:
 * Equilibribum Analysis with Mathematical Programming Methods"
 * Prentice-Hall (out of print, available online from
 * http://web.mit.edu/sheffi/www/urbanTransportation.html
 *   
 ****************************************************************************/

#define TIME_VOLUME_UPDATE
#undef DEBUG_VERIFY_PAPE
#undef TIME_VERIFY_PAPE // mutually exclusive with DEBUG_VERIFY_PAPE

#if defined DEBUG_VERIFY_PAPE && defined TIME_VERIFY_PAPE
#error "cannot have both DEBUG_VERIFY_PAPE and TIME_VERIFY_PAPE"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include "parsetapfiles.h"
#include "pape_cuda_host.h"
#include "tap_functions.h"
#include "utils.h"
#include "volume_update_host.h"
#if defined DEBUG_VERIFY_PAPE || defined TIME_VERIFY_PAPE
#include "sssp_pape.h"
#endif



/*
   In this program, all data for links will be in (1d) arrays such
   that the entries are sorted by link origin (ascending) and then by
   link destinatino (ascending). This way we can have net.links[k]
   being the data about link k and link_volumes[k] and link_costs[k]
   the current volume and cost for link k, the latter can be used
   directly as the weight Wa[k] in the packed adjacency list
   representation of the network (used for efficient shortest paths
   computation) since that is sorted this way.

   Note that node (and zone) numbers start at 1
   (since nodes are numbered in input and output files this way) but
   links start normal C style at the 0th entry in the array as they
   are not numbered externally.  Hence we are always adding 1 to
   num_nodes in allocating vectors.
*/
   
   
/*****************************************************************************
 *
 * constants
 *
 ****************************************************************************/


/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/

/*
 * host_link_volume_update - update the link volume vector (host version)
 *
 * Using the results from the shortest path computatinos, update the volume
 * on each path.
 * This version runs (single-threaded) on host CPU, no atomic operations needed.
 * It is faster to use the multithreaded versino on GPU, 
 * which requires atomic instructions.
 *
 * Parameters:
 *   links - net link data 
 *    num_links - number of links (edges)
 *   num_start_nodes - number of start nodes (=zones)
 *    demands  - demands array parsed from trips file
 *                for each origin demands[origin] is an array of demands structs
 *                terminated by one with 0  dest (origins and dests start at 1)
 *    link_volumes - (OUT) - volume on each link
 *    predlink  -  predecessor link array,
 *               must have space for num_start_nodes*num_nodes entries
 *
 *  Return value;
 *    None.
 *
 */

static void host_link_volume_update(link_data_t links[],
                                    long num_links,
                                    long num_start_nodes,
                                    demand_data_t *demands[],
                                    double link_volumes[],
                                    long predlink[])
{
  long orig,dest;
  long v,k,i;
  double routeflow;

  for (k = 0; k < num_links; k++)
    link_volumes[k] = 0;

  for (orig = 1; orig < num_start_nodes+1; orig++)
  {
//    fprintf(stderr, "orig = %d\n", orig); //XXX

    for (i = 0; (dest = demands[orig][i].dest) != 0; i++)
    {
//      fprintf(stderr, "orig = %d dest = %d\n", orig, dest); // XXX
      routeflow = demands[orig][i].demand;
      if (orig == dest || routeflow <= 0)/*should never be -ve; stops warning*/
        continue;
      v = dest;
      while (v != orig)
      {
        k = predlink[v * num_start_nodes + orig];
        if (k < 0)
        {
          /* this shouldn't happen, it means there is a nonzero demand
             between two nodes but no path between them  - 
             probably due to bad input data */
          fprintf(stderr, "ERROR: no path from %ld to %ld\n", orig, dest);
          exit(1);
          break;
        }
        assert(k >= 0);
        assert(links[k].term_node == v);
        link_volumes[k] += routeflow;
        v = links[k].init_node;
      }
    }
  }
}


/*****************************************************************************
 *
 *  functions
 *
 ****************************************************************************/

/*
 *  fw_assign() - flow assignment in Frank-Wolfe
 *
 *   Assign flow from O-D demands to links according to current shortest
 *   paths (lowest costs) on netgraph
 *
 * Parameters:
 *    net - net structure parsed from net file
 *    demands  - demands array parsed from trips file
 *                for each origin demands[origin] is an array of demands structs
 *                terminated by one with 0  dest (origins and dests start at 1)
 *    Va, Ea, Wa - network in packed adjancey list format (Wa is current costs)
 *    link_volumes - (IN/OUT) - volume on each link
 *   dist (WORK) - vector of num_nodes doubles for distance from origin of each
 *                 for sssp_pape()
 *    predlink  -  (WORK) predecessor link array,
 *               must have space for num_start_nodes*num_nodes entries
 *
 *  Return value;
 *    None.
 */
void fw_assign(net_data_t *net, demand_data_t *demands[],
                      long Va[], long Ea[], double Wa[], 
                      double link_volumes[],
                      double dist[], long predlink[])
{
  long k;
  long *start_nodes = NULL;;
  long num_start_nodes = net->num_zones;

  // TODO  initialize start_nodes only once at start
  // TODO also set up Va, Ea only once at start, only WA changes,but < 1 ms to copy so not much point


  if (!(start_nodes = (long *)malloc(num_start_nodes * sizeof(long))))
  {
    fprintf(stderr, "malloc start_nodes failed\n");
    exit(1);
  }
  for (k = 0; k < num_start_nodes; k++)
    start_nodes[k] = k+1;
 

  pape_cuda(Va, Ea, Wa, net->num_nodes+1, net->num_links, 
            start_nodes, num_start_nodes, net->first_thru_node,
            dist, predlink);

#if defined DEBUG_VERIFY_PAPE || defined TIME_VERIFY_PAPE
  fprintf(stderr, "DEBUG: running sssp_pape() on host to verify...");
  double *gold_dist = (double *)malloc((net->num_nodes+1)*sizeof(double));
  long  *gold_predlink = (long*)malloc((net->num_nodes+1)*sizeof(long));
  long * queue_next = (long*)malloc((net->num_nodes+1)*sizeof(long));
#if defined TIME_VERIFY_PAPE
  struct timeval pape_start_timeval, pape_end_timeval, pape_elapsed_timeval;
  int pape_etime;
  gettimeofday(&pape_start_timeval, NULL);
#elif defined DEBUG_VERIFY_PAPE
  bool failed = false;
#endif
  for (long orig = 1; orig  <=  num_start_nodes; orig++)
  {
    sssp_pape(Va, Ea, Wa, net->num_nodes+1, net->num_links, 
              orig,  net->first_thru_node,
              gold_predlink, gold_dist, queue_next);
#if defined DEBUG_VERIFY_PAPE
    long dest;
//    for (long i = 0; (dest = demands[orig][i].dest) != 0; i++)
//    for (dest = 1; dest <= num_start_nodes; dest++)
    for (dest = 1; dest <= net->num_nodes; dest++)
    {
      if (fabs(gold_dist[dest]- dist[dest * num_start_nodes + orig]) > EPS)
      {
        fprintf(stderr, "ERROR orig %ld dest %ld dist=%f gold_dist=%f\n",
                orig, dest, dist[dest*num_start_nodes+orig], gold_dist[dest]);
        failed = true;
      }
      
      if (gold_predlink[dest] != predlink[dest * num_start_nodes+orig])
      {
        fprintf(stderr, "ERROR orig %ld dest %ld prelink=%ld gold_predlink=%ld\n",
                orig, dest, predlink[dest*num_start_nodes+orig], gold_predlink[dest]);
        failed = true;
      }
    }
#endif /* DEBUG_VERIFY_PAPE */
  }
#if defined DEBUG_VERIFY_PAPE
  if (failed)
    fprintf(stderr, "test FAILED\n");
  else
    fprintf(stderr, "test OK\n");
#elif defined TIME_VERIFY_PAPE
  gettimeofday(&pape_end_timeval, NULL);
  timeval_subtract(&pape_elapsed_timeval, &pape_end_timeval, &pape_start_timeval);
  pape_etime = 1000*pape_elapsed_timeval.tv_sec + pape_elapsed_timeval.tv_usec/1000;
  fprintf(stderr, "sssp_pape host time %d ms\n", pape_etime);
#endif

#endif /* DEBUG_VERIFY_PAPE || TIME_VERIFY_PAPE */

#ifdef TIME_VOLUME_UPDATE
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int etime;
  gettimeofday(&start_timeval, NULL);
#endif

#ifdef USE_CUDA_VOLUME_UPDATE
  link_volume_update(num_start_nodes,net->num_links, net->num_nodes, link_volumes, predlink);
#else
  host_link_volume_update(net->links, net->num_links, num_start_nodes, demands,  link_volumes, predlink);
#endif

#ifdef TIME_VOLUME_UPDATE
  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  etime = 1000*elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
  fprintf(stderr, "total link volume update time %d ms\n", etime);
#endif

  free(start_nodes);
}

