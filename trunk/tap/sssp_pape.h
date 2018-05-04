#ifndef SSSP_PAPE_H
#define SSSP_PAPE_H
/*****************************************************************************
 * 
 * File:    sssp_pape.h
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: sssp_pape.h 699 2011-09-14 04:41:40Z astivala $
 *
 * single-source shortest paths using 
 * the d'Esopo-Pape algorithm [Pape 1974 "Implementation and Efficiency
 * of Moore-Algorithms for the shortest route problem" Math. Prorgam. 7:212-222].
 * This algorithm is well suited to shorest paths in road networks
 * (see Klunder & Post 2006 "The Shortest Path Problem on Large-Scale
 * Real-Road Networks" Networks 48(4):182-194).
 *
 * The d'Esopo-Pape is termed a "label correcting" rather than "label setting"
 * (Dijkstra type) algorithm. A node is always removed from head of queue,
 * and is placed at the tail of the queue if it has never been in the queue
 * before, otehrwised placed at the head.
 *
 ****************************************************************************/

#include <values.h> /* FLT_MAX, DBL_MAX */
#ifdef SOLARIS
#include <float.h>
#endif

#define FLOATINF      DBL_MAX

#ifdef __cplusplus
extern "C" {
#endif

void sssp_pape(long Va[],long Ea[], double Wa[], long num_nodes, long num_edges,
               long start_node, long first_thru_node, long prev[],
               double dist[], long queue_next[]);

void sssp_slf(long Va[],long Ea[], double Wa[],
                   long num_nodes, long num_edges,
                   long start_node, long first_thru_node, long prev[],
                   double dist[], long queue_next[]);

void sssp_slf_lll(long Va[],long Ea[], double Wa[],
                   long num_nodes, long num_edges,
                   long start_node, long first_thru_node, long prev[],
                   double dist[], long queue_next[]);

void sssp_pape_lll(long Va[],long Ea[], double Wa[],
                   long num_nodes, long num_edges,
                   long start_node, long first_thru_node, long prev[],
                   double dist[], long queue_next[]);

void get_shortest_path(long pred[], long orig, long dest, long path[], 
                       long *pathlen);


void sssp_prevnodefirst(long Va[],long Ea[], double Wa[],
                        long num_nodes, long num_edges,
                        long start_node, long first_thru_node, 
                        long old_prevnode[],
                        long prevlink[], long prevnode[],
                        double dist[], long queue_next[]);

#ifdef __cplusplus
}
#endif

#endif /* SSSP_PAPE_H */
