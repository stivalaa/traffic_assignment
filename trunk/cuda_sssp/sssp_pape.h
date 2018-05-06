#ifndef SSSP_PAPE_H
#define SSSP_PAPE_H
/*****************************************************************************
 * 
 * File:    sssp_pape.h
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: sssp_pape.h 260 2011-04-28 04:29:24Z astivala $
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

#define FLOATINF      FLT_MAX


#ifdef __cplusplus
extern "C" {
#endif

void sssp_pape(int Va[],int Ea[], float Wa[], int num_nodes, int num_edges,
               int start_node, int first_thru_node, 
               float distances[], int prev[]);

void sssp_slf(int Va[],int Ea[], float Wa[], int num_nodes, int num_edges,
               int start_node, int first_thru_node, 
               float distances[], int prev[]);

void sssp_pape_lll(int Va[],int Ea[], float Wa[], 
                   int num_nodes, int num_edges,
                   int start_node, int first_thru_node, 
                   float distances[], int prev[]);
  
void sssp_bellmanford(int Va[],int Ea[], float Wa[], int num_nodes, 
              int num_edges,
               int start_node, int first_thru_node, 
               float distances[], int prev[]);

void sssp_slf_lll(int Va[],int Ea[], float Wa[], int num_nodes, int num_edges,
               int start_node, int first_thru_node, 
               float distances[], int prev[]);

void get_shortest_path(int pred[], int orig, int dest, int path[], 
                       int *pathlen);

#ifdef __cplusplus
}
#endif

#endif /* SSSP_PAPE_H */
