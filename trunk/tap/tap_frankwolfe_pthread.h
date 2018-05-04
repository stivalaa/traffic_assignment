#ifndef TAP_FRANKWOLFE_PTHREAD_H
#define TAP_FRANKWOLFE_PTHREAD_H
/*****************************************************************************
 * 
 * File:    tap_frankwolfe_pthread.c
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_frankwolfe_pthread.h 787 2011-10-05 23:41:39Z astivala $
 *   
 *
 ****************************************************************************/

#include "parsetapfiles.h"

/* constants */

#define MAX_NUM_THREADS  512 /* max number of threads we can allow */


/* global data */

extern unsigned int num_cores;   /* number of cores found on system */
extern unsigned int num_threads; /* maximum number of worker threads allowed */


/* function prototypes */


int  tap_frankwolfe(net_data_t *net,
                    demand_data_t **demands,
                    int warm_start_mode,
                    FILE *flows_input_fp,
                    int target_iterations,
                    double target_relgap, FILE *flow_fp,
                    double *total_cost_at_end, double *objvalue_at_end,
		    int *iterations_at_end, double *relgap_at_end,
            int time_each_iteration);


#endif /* TAP_FRANKWOLFE_PTHREAD_H */
