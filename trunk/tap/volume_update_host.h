#ifndef VOLUME_UPDATE_HOST_H
#define VOLUME_UPDATE_HOST_H
/*****************************************************************************
 * 
 * File:    volume_update_host.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: volume_update_host.h 701 2011-09-14 06:24:14Z astivala $
 *
 * CUDA host code for volume vector update
 *
 ****************************************************************************/

#undef USE_CUDA_VOLUME_UPDATE // FIXME this only really works on Fermi and even then will unpredicably cause kernel launch failure


void link_volume_update(long num_start_nodes, long num_edges, long num_nodes,
                        double link_volumes[],
                        long predlink[]);


void link_volume_data_setup(link_data_t *links, long num_links,
                            long num_start_nodes,
                            demand_data_t *demands[]);

void link_volume_data_cleanup(void);


#endif /* VOLUME_UPDATE_HOST_H */
