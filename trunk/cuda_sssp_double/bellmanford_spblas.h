#ifndef BELLMANFORD_SPBLAS_H
#define BELLMANFORD_SPBLAS_H
/*****************************************************************************
 * 
 * File:    bellmanford_spblas.h
 * Author:  Alex Stivala
 * Created: November 2011
 *
 * $Id: bellmanford_spblas.h 835 2011-11-08 05:46:16Z astivala $
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
 ****************************************************************************/

#include <values.h> /* FLT_MAX, DBL_MAX */
#include "sssp.h"
#define FLOATINF      DBL_MAX

#ifdef __cplusplus
extern "C" {
#endif


void bellmanford_spblas(adjlist_entry_t *adjlist, 
                        long num_nodes, 
                        long num_edges,
                        long start_node,
                        double dist[]);

#ifdef __cplusplus
}
#endif

#endif /* BELLMANFORD_SPBLAS_H */
