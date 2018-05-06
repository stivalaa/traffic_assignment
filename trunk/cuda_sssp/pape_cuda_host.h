#ifndef PAPE_CUDA_HOST_H
#define PAPE_CUDA_HOST_H
/*****************************************************************************
 * 
 * File:    pape_cuda_host.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: pape_cuda_host.h 221 2011-04-12 07:13:04Z astivala $
 *
 * CUDA host code for CUDA implemnetatnion of d'Esopo-Pape algorithm.
 *
 ****************************************************************************/

void pape_cuda(int Va[], int Ea[], float Wa[], 
               int num_nodes, int num_edges,
               int start_nodes[], int num_start_nodes,
               float *distances, int *predecessors);

#endif /* PAPE_CUDA_HOST_H */
