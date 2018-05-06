#ifndef SSSP_GOLD_H
#define SSSP_GOLD_H
/*****************************************************************************
 * 
 * File:    sssp_gold.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: sssp_gold.h 668 2011-09-08 04:40:08Z astivala $
 *
 ****************************************************************************/

#include "sssp.h"

double sssp_gold(adjlist_entry_t adjlist[], long num_nodes, long num_arcs,
               long source,
               double distances[], long predecessors[]);

#endif /* SSSP_GOLD_H */
