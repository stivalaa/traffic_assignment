/*****************************************************************************
 * 
 * File:    bellmanford_spblas.c
 * Author:  Alex Stivala
 * Created: November 2011
 *
 * $Id: bellmanford_spblas.c 849 2011-11-14 21:42:09Z astivala $
 *
 * single-source shortest paths using  algebraic Bellman-Ford algorithm
 * implemented with Sparse BLAS (modified to use (min, +) algebra rather
 * than usual (+, *))
 *
 * See:
 * Fineman and Robinsion "Fundamental Graph Algorithmns", Chapter 5 in
 * Kepner and Gilbert (eds) "Graph Algorithms in the Language of Linear Algebra"
 * SIAM, 2011, ISBN 978-0-898719-90-1
 *
 *
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

#include "blas_sparse.h"

#include "bellmanford_spblas.h"


#ifdef DEBUG
static void dumpvec(double v[], int n)
{
  int i;
  for (i = 0; i < n; i++)
    printf("%4.2f ", v[i]);
  printf("\n");
    }
#endif /* DEBUG */


/*
 * bellmanford_spblas() - single-source shortest path by algebraic Bellman-Ford 
 *  
 *
 * This is the algebraic formulation of Bellmnan-Ford algorithm,  with loop over
 * all nodes, computing matrix-vector multiplication in (min, +) algebra,
 * using modified Sparse BLAS reference implemntation.
 * 
 *
 * Parameters:
 *    adjlist - graph in adjacency list format
 *    num_nodes - number of nodes 
 *    num_edges - number of edges (elements in adjlist)
 *    start_node - source node
 *    dist - (OUT) distance array, dist from source to each node
 */

void bellmanford_spblas(adjlist_entry_t *adjlist, 
                        long num_nodes, 
                        long num_edges,
                        long start_node,
                        double d[])
{
  long i;
  double alpha = 0; /* for (min,+) algebra this is the 'multiplicative' (+) identity */
  int n = (int)num_nodes;
  blas_sparse_matrix A;
  double *y;

  /* allocate space for y vector */
  if (!(y = (double *)malloc(n * sizeof(double))))
  {
    fprintf(stderr, "malloc vector failed\n");
    exit(1);
  }
  
  /* build the adjacency sparse matrix */
  A = BLAS_duscr_begin(n, n);

  for (i = 0; i <n; i++)
    BLAS_duscr_insert_entry(A, 0, i, i); /* 'zero' ('additive' identity) in this algebra is inf */
                                         /* but we want 'multiplicative' identity on diagonal */

  for (i = 0; i < num_edges; i++)
  {
/*    fprintf(stderr, "XXX n = %d i = %d from = %d to = %d\n", n,i, (int)adjlist[i].from,(int)adjlist[i].to); */
    BLAS_duscr_insert_entry(A, adjlist[i].cost, (int)adjlist[i].from, (int)adjlist[i].to);
  }
  BLAS_duscr_end(A);

  /* initialize dist vector, all inf except start node 0 */
  /* the output accumulator dist vector y must be all 'additive' identity inf  (for min)*/
  for (i = 0; i < num_nodes; i++)
  {
    d[i] = FLOATINF;
    y[i] = FLOATINF;
  }
  d[start_node] = 0;

  

  /* compute d = (A^T)*d i.e. d = d*A
     for each node where * is matrix-vector mult in (min,+) algebra */
  for (i = 0; i < num_nodes; i++)
  {
    BLAS_dusmv(blas_trans, alpha, A, d, 1, y, 1);
    memcpy(d, y, num_nodes*sizeof(double));
#ifdef DEBUG
/*    dumpvec(d, n); */
#endif

  }

  /* cleanup Sparse BLAS  matrix handle and free memory */
  BLAS_usds(A);
  free(y);
}
