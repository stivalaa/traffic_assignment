#ifndef HARISH_HOST_H
#define HARISH_HOST_H
/*****************************************************************************
 * 
 * File:    harish_host.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: harish_host.h 97 2011-02-21 06:51:13Z astivala $
 *
 * CUDA host code for single-source shortest path implementation using
 * CUDA based on:
 *
 * Harish and Narayanan 2007 "Accelerating Large Graph Algorithms on the GPU
 * Using CUDA" HiPC 2007, LNCS 4873: 197-208
 * 
 *
 ****************************************************************************/

void harish_sssp(int Va[], int Ea[], float Wa[], int alloc_num_nodes,
                 int num_nodes, int num_edges, 
                 int start_node, float distances[], int predecessors[]);

#endif /* HARISH_HOST_H */
