#ifndef VOLUME_UPDATE_KERNELS_H
#define VOLUME_UPDATE_KERNELS_H
/*****************************************************************************
 * 
 * File:    volume_update_kernels.h
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: volume_update_kernels.h 692 2011-09-13 07:18:19Z astivala $
 *
 *
 *
 ****************************************************************************/

#include "tap_types.h" /* only for link_data_t */

__global__ void init_volume_kernel(double link_volumes[], long num_edges);

__global__ void link_volume_update_kernel(long link_init_nodes[],
                                          long num_start_nodes,
                                          double od_matrix[],
                                          double link_volumes[],
                                          long predlink[]);
                                          

#endif
