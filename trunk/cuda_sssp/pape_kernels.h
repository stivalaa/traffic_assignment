#ifndef PAPE_KERNELS_H
#define PAPE_KERNELS_H
/*****************************************************************************
 * 
 * File:    pape_kernels.h
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: pape_kernels.h 221 2011-04-12 07:13:04Z astivala $
 *
 * CUDA implementation of single-source shortest paths using the
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
 ****************************************************************************/

__global__ void pape_init_arrays(int num_nodes, int num_start_nodes,
                                 int start_nodes[],
                                 float dist[], int prev[], int queue_next[]);

__global__ void pape_kernel(int Va[],int Ea[], float Wa[],
                            int num_nodes, int num_edges,
                            int num_start_nodes,
                            int start_nodes[], 
                            int first_thru_node, float dist[], int prev[],
                            int queue_next[]);


#endif /* PAPE_KERNELS_H */
