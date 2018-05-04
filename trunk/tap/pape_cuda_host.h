#ifndef PAPE_CUDA_HOST_H
#define PAPE_CUDA_HOST_H
/*****************************************************************************
 * 
 * File:    pape_cuda_host.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: pape_cuda_host.h 690 2011-09-13 06:42:49Z astivala $
 *
 * CUDA host code for CUDA implemnetatnion of d'Esopo-Pape algorithm.
 *
 ****************************************************************************/

void pape_cuda(long Va[], long Ea[], double Wa[], 
               long num_nodes, long num_edges,
               long start_nodes[], long num_start_nodes,
               long first_thru_node,
               double *distances, long *predecessors);

#endif /* PAPE_CUDA_HOST_H */
