#ifndef TAP_FUNCTIONS_H
#define TAP_FUNCTIONS_H
/*****************************************************************************
 * 
 * File:    tap_functions.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_functions.h 671 2011-09-08 09:03:36Z astivala $
 *
 * Traffic assignment functions
 *
 *
 ****************************************************************************/

#include <values.h> /* FLT_MAX, DBL_MAX */
#ifdef SOLARIS
#include <float.h>
#endif

#define FLOATINF      DBL_MAX

static const double EPS = 1e-08; /* epsilon for difference between fp numbers */
static const int NUM_STEPS = 10; /* number of steps in cost function */

#ifdef __cplusplus
extern "C" {
#endif

 double bpr_cost_function(const link_data_t *link, double Q);
 double linkcost_integral(const link_data_t *link, double Q);
 double links_objective_function(const net_data_t *net, double volumes[]);
 double link_directional_derivative(const net_data_t *net,
                                     double volume[],
                                     double delta_volume[], double lambda);
 double cost_step(link_data_t *link, double vol);
 double cost_step_function(link_data_t *link, double vol);
 double distance_from_next_step(link_data_t links[], 
                                      double link_volumes[],
                               long link_idx);
 int adjlist_entry_compar(const void *ent1, const void *ent2);
 void adjlist_to_packed_arrays(link_data_t adjlist[],
                              long num_nodes,
                              long num_edges,
                              long Va[],
                              long Ea[],
                              double Wa[]);

 void update_link_costs(net_data_t *net, double link_volumes[],
                       double link_costs[]);
 void update_link_step_costs(net_data_t *net, double link_volumes[],
                       double link_costs[]);

 double demand_sum(demand_data_t * const demands[], long num_zones);
 double total_link_cost(const net_data_t *net,
                              const double link_volumes[],
                       const double link_costs[]);
 void print_flow(FILE *fp, const net_data_t *net, 
                       const double link_volumes[],
                const double link_costs[]);

 double line_search(const net_data_t *net,
                    double link_volumes[], double link_delta_volumes[]);

 int apply_trip_modification(demand_data_t *demands[], trip_mod_t *mod,
                              long num_zones);

#ifdef __cplusplus
}
#endif

#endif /* TAP_FUNCTIONS_H */
