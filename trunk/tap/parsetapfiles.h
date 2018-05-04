#ifndef PARSETAPFILES_H
#define PARSETAPFILES_H
/*****************************************************************************
 * 
 * File:    parsetapfiles.h
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: parsetapfiles.h 682 2011-09-12 06:40:52Z astivala $
 *
 * Functions to parse the test data _trips, _net, (traffic assignment input)
 * data in the format from
 *
 * http://www.bgu.ac.il/~bargera/tntp/
 *
 ****************************************************************************/


#include "tap_types.h"

/****************************************************************************
 *
 * function prototypes
 *
 ****************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

/* parse network data */
int parse_net_file(FILE *fp, net_data_t *net);

/* parse demand (trips) data */
int parse_trips_file(FILE *fp, demand_data_t **demands[], long *num_zones);

/* parse flows (TAP output) data */
int parse_flows_file(FILE *flow_fp, link_flow_t **link_flows, long *num_links);

/* parse netowkr modifications file */
int parse_net_mod_file(FILE *fp, net_mod_t **mods, int *num_mods);

/* qsort comparison function for link data  in network data */
int link_data_compar(const void *ent1, const void *ent2);

/* parse trip data modifications file */
int parse_trip_mod_file(FILE *fp, trip_mod_t **mods, int *num_mods);

#ifdef __cplusplus
}
#endif


#endif /* PARSETAPEFILES_H */
