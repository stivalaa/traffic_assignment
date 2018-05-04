#ifndef TAP_FRANKWOLFE_CUDA_H
#define TAP_FRANKWOLFE_CUDA_H
/*****************************************************************************
 * 
 * File:    tap_frankwolfe_cuda.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_frankwolfe_cuda.h 672 2011-09-08 09:04:30Z astivala $
 *   
 *
 ****************************************************************************/

#include "parsetapfiles.h"

/* constants */


/* global data */

extern int using_fermi_architecture; // using a Fermi architecture GPU


/* function prototypes */

void fw_assign(net_data_t *net, demand_data_t *demands[],
                      long Va[], long Ea[], double Wa[], 
                      double link_volumes[],
               double dist[], long predlink[]);

#endif /* TAP_FRANKWOLFE_CUDA_H */
