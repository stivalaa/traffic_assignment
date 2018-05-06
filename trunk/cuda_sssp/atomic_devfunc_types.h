#ifndef ATOMIC_DEVFUNC_TYPES_H
#define ATOMIC_DEVFUNC_TYPES_H
/*****************************************************************************
 * 
 * File:    atomic_devfunc_types.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: atomic_devfunc_types.h 209 2011-04-10 23:33:11Z astivala $
 *
 * CUDA atomic operation device functions.
 *
 * Requires device compute capabillity at least 1.2 (uses 64 bit atomic
 * functions).
 *
 ****************************************************************************/

/****************************************************************************
 * 
 * type definitions
 *
 ****************************************************************************/

/* cost and node together so we can operate on them with single
   64 bit atomic instruction */
typedef struct __align__(8) cost_node_pair_s  // must be aligned as long long
{
    float cost;
    int   node;
} cost_node_pair_t;

#endif /* ATOMIC_DEVFUNC_TYPES */
