#ifndef OKUYAMA_HOST_H
#define OKUYAMA_HOST_H
/*****************************************************************************
 * 
 * File:    okuyama_host.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: okuyama_host.h 99 2011-02-21 23:25:38Z astivala $
 *
 *
 * CUDA host code for multiple-source shortest path implementation using
 * CUDA based on:
 *
 * Okuyama, Ino, Hagihara 2008 "A Task Parallel Algorithm for Computing
 * the costs of All-Pairs Shortest Paths on the CUDA-compatible GPU"
 * Intl. Symp. Parallel Distributed Processing with Applications. 284-291
 *
 ****************************************************************************/

void okuyama_sssp(int Va[], int Ea[], float Wa[], 
                 int num_nodes, int num_edges, 
                  int start_nodes[], int num_start_nodes, float distances[],
                  int predecessors[]);

#endif /* OKUYAMA_HOST_H */
