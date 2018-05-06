#ifndef SSSP_GOLD_H
#define SSSP_GOLD_H
/*****************************************************************************
 * 
 * File:    sssp_gold.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: sssp_gold.h 129 2011-03-02 04:20:36Z astivala $
 *
 ****************************************************************************/

#include "sssp.h"

double sssp_gold(adjlist_entry_t adjlist[], int num_nodes, int num_arcs,
               int source,
               float distances[], int predecessors[]);

#endif /* SSSP_GOLD_H */
