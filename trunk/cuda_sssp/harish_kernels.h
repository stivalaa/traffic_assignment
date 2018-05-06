#ifndef HARISH_KERNELS_H
#define HARISH_KERNELS_H
/*****************************************************************************
 * 
 * File:    harish_kernels.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: harish_kernels.h 224 2011-04-13 06:09:10Z astivala $
 *
 * Single-source shortest path implementation using CUDA based on:
 *
 * Harish and Narayanan 2007 "Accelerating Large Graph Algorithms on the GPU
 * Using CUDA" HiPC 2007, LNCS 4873: 197-208
 * 
 *
 * Developed on CUDA 2.3.
 * Requires device compute capabillity at least 1.2 (uses 64 bit atomic
 * functions).
 *
 ****************************************************************************/

//#include "atomic_devfunc.h"

const int BLOCK_SIZE = 128; // number of threads in a block


__global__ void init_mask_cost_update_arrays(bool Ma[], 
                                             float Ca[],
                                             struct cost_node_pair_s Ua[], 
                                             int s,
                                             int   Pa[],
                                             int num_nodes);

__global__ void  harish_scatter_kernel(int Va[],
                                       int Ea[],
                                       float Wa[],
                                       bool Ma[],
                                       float Ca[],
                                       volatile struct cost_node_pair_s Ua[],
                                       int Pa[],
                                       int num_nodes);

__global__ void harish_update_kernel(int Va[], int Ea[], float Wa[],
                                     bool Ma[], float Ca[], 
                                     struct cost_node_pair_s Ua[],
                                     int Pa[], 
                                     int num_nodes);

#endif /* HARISH_KERNELS_H */

