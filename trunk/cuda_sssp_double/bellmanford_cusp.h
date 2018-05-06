#ifndef BELLMANFORD_CUSP_H
#define BELLMANFORD_CUSP_H
/*****************************************************************************
 * 
 * File:    bellmanford_cusp.h
 * Author:  Alex Stivala
 * Created: November 2011
 *
 * $Id: bellmanford_cusp.h 840 2011-11-09 05:50:30Z astivala $
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
 ****************************************************************************/

#include <values.h> /* FLT_MAX, DBL_MAX */
#include "sssp.h"
#define FLOATINF      DBL_MAX

#ifdef __cplusplus
extern "C" {
#endif


void bellmanford_cusp_device(adjlist_entry_t *adjlist, 
                        long num_nodes, 
                        long num_edges,
                        long start_node,
                        double dist[]);

void bellmanford_cusp_host(adjlist_entry_t *adjlist, 
                        long num_nodes, 
                        long num_edges,
                        long start_node,
                        double dist[]);

#ifdef __cplusplus
}
#endif

#endif /* BELLMANFORD_CUSP_H */
