/*****************************************************************************
 * 
 * File:    pape_kernels.c
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: pape_kernels.cu 256 2011-04-27 07:20:26Z astivala $
 *
 * CUDa implementation ofsingle-source shortest paths using the
 * d'Esopo-Pape algorithm [Pape 1974 "Implementation and Efficiency of
 * Moore-Algorithms for the shortest route problem"
 * Math. Prorgam. 7:212-222].  This algorithm is well suited to
 * shorest paths in road networks (see Klunder & Post 2006 "The
 * Shortest Path Problem on Large-Scale Real-Road Networks" Networks
 * 48(4):182-194).
 *
 * The d'Esopo-Pape is termed a "label correcting" rather than "label setting"
 * (Dijkstra type) algorithm. A node is always removed from head of queue,
 * and is placed at the tail of the queue if it has never been in the queue
 * before, otehrwised placed at the head.
 *
 * TODO SLF and LLL variations (as per sssp_pape.c) and bidirectional.
 *
 ****************************************************************************/


#include <cutil_inline.h>      /* CUDA SDK */

#include "pape_kernels.h"


/****************************************************************************
 *
 * constants and type definitions
 *
 ****************************************************************************/

#define INVALID      -1
#define WAS_IN_QUEUE -2


/****************************************************************************
 *
 * __global__ functions: GPU kernels, callable from host
 *
 ****************************************************************************/


/*
 * pape_init_arrays() - kernel to init arrays for pape_kernel()
 *
 * Parameters:
 *    num_nodes - number of nodes (elemnts in  dist, queue_next)
 *    num_start_nodes - number of source nodes
 *    start_nodes - source nodes
 *    dist - (OUT) distance array, dist from each source to each node
 *                  dist[i,s] is dist from s to i
 *                   must be intiizlize dto all INVALID except start_node=0
 *    prev  -  (OUT) predecessor  array,
 *               must have space for num_nodes*num_setart_nodes entries
 *              Each prev[i,s] is predecessor of node i for startnode s
 *              in shortest path to node i from s
 *                intiizlied to all INVALID
 *    queue_next (OUT]) - array of length num_nodes for queue
 *                         initilized to all INVALID
 */
__global__ void pape_init_arrays(int num_nodes, int num_start_nodes,
                                 int start_nodes[],
                                 float dist[], int prev[], int queue_next[])
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
      prev[v * num_start_nodes + s] = INVALID;
      queue_next[v * num_start_nodes + s] = INVALID;
      if (v == s)
        dist[v * num_start_nodes + s] = 0;
      else
        dist[v * num_start_nodes + s] = FLOATINF;
    }
  }
}


/*
 * pape_kernel() - single-source shortest path by d'Esopo-Pape algorithm
 *
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    num_start_nodes - number of source nodes
 *    start_node - source nodes
 *    first_thru_node - first node number that is allowed in a path
 *                      (earlier ones are actually 'zones' for origin/dest).
 *    dist - (in/OUT) distance array, dist from source to each node
 *                   dist[s,i] is dist from s to i
 *                   must be intiizlize dto all INVALID except start_node=0
 *    prev  -  (OUT) predecessor  array,
 *               must have space for num_nodes entries
 *              Each prev[i,s] is predecessor of node i 
 *              in shortest path to node i from source s
 *               must be intiizlied to all INVALID
 *    queue_next ([workspace]) - array of length num_nodes for queue
 *                        must be initilized to all INVALID
 */
__global__ void pape_kernel(int Va[],int Ea[], float Wa[],
                            int num_nodes, int num_edges,
                            int num_start_nodes,
                            int start_nodes[], 
                            int first_thru_node, float dist[], int prev[],
                            int queue_next[])
{

  int u,v;
  int i;
  int queue_first, queue_last;
  float uvdist, newcost;


  // each thread does one start node (not caring about block/grid)
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < num_start_nodes)
  {
  int s = start_nodes[idx];

//  assert (!(start_node >= num_nodes));

  queue_first = INVALID;
  queue_last = INVALID;

  u = s;
  while (u != INVALID && u != WAS_IN_QUEUE)
  {
    if (u >= first_thru_node || u == s)
    {
      for (i = Va[u]; i < Va[u+1]; i++) /* all neighbours of u */
      {
        v = Ea[i];      /* v is adjacent to u */
//        assert(v >= 0 && v < num_nodes);
        uvdist = Wa[i]; /* with this cost on edge */
        newcost = dist[u * num_start_nodes + s] + uvdist;
        if (newcost < dist[v * num_start_nodes + s])
        {
          dist[v * num_start_nodes + s] = newcost;
          prev[v * num_start_nodes + s] = u;
          if (queue_next[v * num_start_nodes + s] == WAS_IN_QUEUE)
          {
            queue_next[v * num_start_nodes + s] = queue_first;
            queue_first = v;
            if (queue_last == INVALID)
            {
              queue_last = v;
            }
          }
          else if (queue_next[v * num_start_nodes + s] == INVALID &&
                   v != queue_last)
          {
            if (queue_last != INVALID)
            {
              queue_next[queue_last * num_start_nodes + s] = v;
              queue_last = v;
            }
            else
            {
              queue_first = v;
              queue_last = v;
              queue_next[queue_last * num_start_nodes + s] = INVALID;
            }
          }
        }
      }
    }
    u = queue_first;
    if (u == INVALID || u == WAS_IN_QUEUE)
      break;
//    assert(u >=0 && u < num_nodes);
    queue_first = queue_next[u * num_start_nodes + s];
    queue_next[u * num_start_nodes + s] = WAS_IN_QUEUE;
    if (queue_last == u)
      queue_last = INVALID;
  }
  } 
}

