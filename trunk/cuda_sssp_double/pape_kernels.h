#ifndef PAPE_KERNELS_H
#define PAPE_KERNELS_H
/*****************************************************************************
 * 
 * File:    pape_kernels.h
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: pape_kernels.h 668 2011-09-08 04:40:08Z astivala $
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

__global__ void pape_init_arrays(long num_nodes, long num_start_nodes,
                                 long start_nodes[],
                                 double dist[], long prev[], long queue_next[]);

__global__ void pape_kernel(long Va[],long Ea[], double Wa[],
                            long num_nodes, long num_edges,
                            long num_start_nodes,
                            long start_nodes[], 
                            long first_thru_node, double dist[], long prev[],
                            long queue_next[]);


#endif /* PAPE_KERNELS_H */
