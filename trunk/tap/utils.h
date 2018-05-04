#ifndef UTILS_H
#define UTILS_H
/*****************************************************************************
 * 
 * File:    utils.h
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * Miscellaneous utilty functions
 *
 * $Id: utils.h 671 2011-09-08 09:03:36Z astivala $
 *
 ****************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

int iDivUp(int a, int b);
int get_num_cores(void);
int timeval_subtract (struct timeval *result, struct timeval *x, 
                       struct timeval *y);



#ifdef __cplusplus
}
#endif

#endif /* UTILS_H */

