/*****************************************************************************
 * 
 * File:    bellmanford_cusp.cpp
 * Author:  Alex Stivala
 * Created: November 2011
 *
 * $Id: bellmanford_cusp.cu 846 2011-11-10 04:28:52Z astivala $
 *
 * single-source shortest paths using algebraic Bellman-Ford algorithm
 * implemented with CUSP sparse matrix library spmv (modified to use (min,
 * +) algebra rather than usual (+, *))
 *
 * See:
 * Fineman and Robinsion "Fundamental Graph Algorithmns", Chapter 5 in
 * Kepner and Gilbert (eds) "Graph Algorithms in the Language of Linear Algebra"
 * SIAM, 2011, ISBN 978-0-898719-90-1
 *
 *
 * for CUSP see http://code.google.com/p/cusp-library/
 * bu tNB we neeed to modify the CUSP code to use (min,+) algebra not usual
 * (+,*).
 *                                                                      
 *
 * This module is compiled with two versions, if CUDSP_DEVICE macro is defined
 * then bellmanford_cusp_device() is built whic uses CUDA device otherwise
 * bellmanford_cusp_host() is built which uses host not device.
 *
 * TODO: implement path retrieval, at the moment we only get the distances
 *
 ****************************************************************************/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>

#undef CUSP_USE_TEXTURE_MEMORY // CUSP ELL spmv faster WITHOUT texture memory e.g.  rmat.txt chicagoregionalflows 682.508972 ms without texture, 823.638000 with texture (for 1 start node)

#include "cusp/coo_matrix.h"
//#include "cusp/dia_matrix.h"
#include "cusp/ell_matrix.h"
#include "cusp/multiply.h"
#include "cusp/transpose.h"
#ifdef DEBUG
#include "cusp/print.h"
#endif

#include "bellmanford_cusp.h"




/*
 * bellmanford_cusp() - single-source shortest path by algebraic Bellman-Ford 
 *  
 *
 * This is the algebraic formulation of Bellmnan-Ford algorithm,  with loop over
 * all nodes, computing matrix-vector multiplication in (min, +) algebra,
 * using modified CUSP spmv 
 * 
 *
 * Parameters:
 *    adjlist - graph in adjacency list format
 *    num_nodes - number of nodes 
 *    num_edges - number of edges (elements in adjlist)
 *    start_node - source node
 *    dist - (OUT) distance array, dist from source to each node
 */

#ifdef CUSP_DEVICE
void bellmanford_cusp_device
#else
void bellmanford_cusp_host
#endif

(adjlist_entry_t *adjlist, 
                        long num_nodes, 
                        long num_edges,
                        long start_node,
                      double dist[])
{
  long i;
  int n = (int)num_nodes;

#ifdef CUSP_DEVICE
#define CUSP_MEMORY cusp::device_memory
#else
#define CUSP_MEMORY cusp::host_memory
#endif

  /* build the adjacency sparse matrix in COO format on host and
     copy to device in DIA format (faster in general for spmv) */
  cusp::coo_matrix<int, double, cusp::host_memory>A_coo(n, n, num_edges+n);
  for (i = 0; i <n; i++)
  {
    /* 'zero' ('additive' identity) in this algebra is inf */
    /* but we want 'multiplicative' identity on diagonal */
    A_coo.row_indices[i] = i;
    A_coo.column_indices[i] = i;
    A_coo.values[i] = 0;
  }
  for (i = 0; i < num_edges; i++)
  {
    A_coo.row_indices[n+i] = adjlist[i].from;
    A_coo.column_indices[n+i] = adjlist[i].to;
    A_coo.values[n+i] = adjlist[i].cost;
  }
  /* copy to device in DIA format */
  //cusp::dia_matrix<int, double,cusp::device_memory> A = A_coo; 
   /* use ELL format not DIA as DIA  gets 'what():  dia_matrix fill-in would exceed maximum tolerance' error on larger data */
  cusp::ell_matrix<int, double,CUSP_MEMORY> A = A_coo; 

  // get A transpose
  cusp::ell_matrix<int, double,CUSP_MEMORY> At;
  cusp::transpose(A, At);

  /* initialize dist vector, all inf except start node 0 */
  cusp::array1d<double, CUSP_MEMORY>y(n);
  cusp::array1d<double, CUSP_MEMORY>d(n);
  thrust::fill(d.begin(), d.end(), double(FLOATINF));
  d[start_node] = 0;
 

#ifdef DEBUG
   cusp::print(At);
   cusp::print(d); 
#endif

  /* compute d = (A^T)*d i.e. d = d*A
     for each node where * is matrix-vector mult in (min,+) algebra */
  for (i = 0; i < num_nodes; i++)
  {
    cusp::multiply(At, d, y); /* y <- A*d */
    d = y;
  }

  for (i = 0; i < n; i++)
    dist[i] = d[i];

}
