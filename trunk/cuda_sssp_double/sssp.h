#ifndef SSSP_H
#define SSSP_H
/*****************************************************************************
 * 
 * File:    sssp.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: sssp.h 668 2011-09-08 04:40:08Z astivala $
 *
 ****************************************************************************/

#include <values.h> // FLT_MAX, DBL_MAX

#define FLOATINF DBL_MAX


/* adjacency list entry */
typedef struct adjlist_entry_s 
{
    long from;          /* from node id (0..n-1) */
    long to;            /* to node id (0..n-1) */
    double cost;        /* cost of arc from 'from' to 'to' >= 0 */
} adjlist_entry_t;

#endif /* SSSP_H */

