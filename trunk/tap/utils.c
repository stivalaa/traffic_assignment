/*****************************************************************************
 * 
 * File:    utils.c
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * Miscellaneous utilty functions
 *
 * $Id: utils.c 277 2011-05-10 07:44:10Z astivala $
 *
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include "utils.h"

/*
 * iDivUp(a,b) = ceil(a / b) 
 */
int iDivUp(int a, int b)
{
  return ((a % b) != 0) ? (a / b + 1) : (a / b);
}


/*
 * Return the number of processors online
 */
int get_num_cores(void)
{
  long ncores;
  if ((ncores = sysconf(_SC_NPROCESSORS_ONLN)) < 0)
  {
    fprintf(stderr,  "sysconf() failed (%d)\n", errno);
    exit(1);
  }
  return (int)ncores;
}


/* Subtract the `struct timeval' values X and Y,
   storing the result in RESULT.
   Return 1 if the difference is negative, otherwise 0.  
(from GNU libc manual) */
     
int
timeval_subtract (struct timeval *result, struct timeval *x, 
                  struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }
     
  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;
     
  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

