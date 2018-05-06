#ifndef SSSP_H
#define SSSP_H
/*****************************************************************************
 * 
 * File:    sssp.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: sssp.h 91 2011-02-21 04:45:26Z astivala $
 *
 ****************************************************************************/

#include <values.h> // FLT_MAX

#define FLOATINF FLT_MAX


/* adjacency list entry */
typedef struct adjlist_entry_s 
{
    int from;          /* from node id (0..n-1) */
    int to;            /* to node id (0..n-1) */
    float cost;        /* cost of arc from 'from' to 'to' >= 0 */
} adjlist_entry_t;

#endif /* SSSP_H */

