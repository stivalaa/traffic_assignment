#ifndef OKUYAMA_KERNELS_H
#define OKUYAMA_KERNELS_H
/*****************************************************************************
 * 
 * File:    okuyama_kernels.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: okuyama_kernels.h 224 2011-04-13 06:09:10Z astivala $
 *
 *
 * CUDA kernels for multiple-source shortest path implementation using
 * CUDA based on:
 *
 * Okuyama, Ino, Hagihara 2008 "A Task Parallel Algorithm for Computing
 * the costs of All-Pairs Shortest Paths on the CUDA-compatible GPU"
 * Intl. Symp. Parallel Distributed Processing with Applications. 284-291
 * 
 * Developed on CUDA 2.3.
 * Requires device compute capabillity at least 1.2 (uses 64 bit atomic
 * functions).
 *
 ****************************************************************************/

#include "atomic_devfunc_types.h"

#define MAX_CONSTANT_NODES 16384 /* to fit in 64 KB constant memory */

__global__ void okuyama_init_mask_cost_update_arrays(bool Ma[], 
                                                     float Ca[],
                                                     cost_node_pair_t Ua[], 
                                                     int start_nodes[],
                                                     int num_nodes,
                                                     int num_start_nodes,
                                                     int d_Pa[]);

__global__ void  okuyama_scatter_kernel(int Va[],
                                       int Ea[],
                                       float Wa[],
                                       bool Ma[],
                                       float Ca[],
                                       volatile cost_node_pair_t Ua[],
                                        int start_nodes[],
                                        int num_nodes, int num_start_nodes,
                                        int d_Pa[]);

__global__ void okuyama_update_kernel(
                                     bool Ma[], float Ca[], 
                                      cost_node_pair_t Ua[],
                                      int start_nodes[],
                                      int num_nodes, int num_start_nodes,
                                      int d_Pa[]);

#endif /* OKUYAMA_KERNELS_H */

