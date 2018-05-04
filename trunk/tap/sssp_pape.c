/*****************************************************************************
 * 
 * File:    sssp_pape.c
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: sssp_pape.c 479 2011-07-13 07:14:16Z astivala $
 *
 * single-source shortest paths using 
 * the d'Esopo-Pape algorithm [Pape 1974 "Implementation and Efficiency
 * of Moore-Algorithms for the shortest route problem" Math. Prorgam. 7:212-222].
 * This algorithm is well suited to shorest paths in road networks
 * (see Klunder & Post 2006 "The Shortest Path Problem on Large-Scale
 * Real-Road Networks" Networks 48(4):182-194).
 *
 * The d'Esopo-Pape is termed a "label correcting" rather than "label setting"
 * (Dijkstra type) algorithm. A node is always removed from head of queue,
 * and is placed at the tail of the queue if it has never been in the queue
 * before, otehrwised placed at the head.
 *
 *
 * TODO add bidirectional search, should almost 1/2 time (see Klunder & Post '06)
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "sssp_pape.h"
#include <math.h>

#undef DEBUG
#undef TIME_DEBUG
#undef ITER_DEBUG
#undef SLOW_INTERNAL_DEBUG

#undef USE_KAHAN_SUM /* Use the Kahan summation algorithm for running total */
                      /* total_cost in sssp_pape_lll() */


/****************************************************************************
 *
 * Constants and type definitions
 *
 ****************************************************************************/

#define INVALID      -1
#define WAS_IN_QUEUE -2


/****************************************************************************
 *
 * Local Functions
 *
 ****************************************************************************/


#ifdef DEBUG
/*
 * dump_packed_arrays() -dump packed array representation for debugging
 *
 * Parameters:
 *
 *        Va[num_nodes] -array of indices to head of each adj list in Ea 
 *        Ea[num_edges] -each entry is 'to' node in list for 'from' node
 *        Wa[num_edges] -each entry is cost of corresponding Ea entry
 *        num_nodes  -number of nodes
 *        num_edges  -number of edges
 *
 * Return value:
 *        None.
 */
static void dump_packed_arrays(const long Va[], const long Ea[],
                               const double Wa[],
                               long num_nodes, long num_edges)
{
  int i;
  for (i = 0; i < num_nodes; i++)
    printf("%ld ", Va[i]);
  printf("\n");
  for (i = 0; i < num_edges; i++)      
    printf("%4ld ", Ea[i]);
  printf("\n");
  for (i = 0; i < num_edges; i++)
    printf("%4.2f ", Wa[i]);
  printf("\n");
}
#endif  /* DEBUG */

/****************************************************************************
 *
 * single source shortest path function 
 *
 ****************************************************************************/

/*
 * sssp_pape() - single-source shortest path by d'Esopo-Pape algorithm
 *
 * Code based on the implemenation
 * from  Bar-Gera's sample Frank-Wolfe implemetatnion at
 * http://www.bgu.ac.il/~bargera/tntp/
 *
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    start_node - source node
 *    first_thru_node - first node number that is allowed in a path
 *                      (earlier ones are actually 'zones' for origin/dest).
 *    prevlink  -  (OUT) predecessor link array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is link index of edge into node i
 *              in shortest path to node i
 *   dist (WORK) - vector of num_nodes doubles for distance from origin of each
 *   queue_next (WORK) - vector of num_nodes longs for queue of nodes
 */

void sssp_pape(long Va[],long Ea[], double Wa[], long num_nodes, long num_edges,
               long start_node, long first_thru_node, long prevlink[],
               double dist[], long queue_next[])
{

  long u,v;
  long i;
  long queue_first, queue_last;
  double uvdist, newcost;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif
#ifdef ITER_DEBUG
  int itercount =0;
#endif


  assert (!(start_node >= num_nodes));
  
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prevlink[v] = INVALID;
    queue_next[v] = INVALID;
    if (v == start_node)
      dist[v] = 0;
    else
      dist[v] = FLOATINF;
  }
  queue_first = INVALID;
  queue_last = INVALID;

  u = start_node;
  while (u != INVALID && u != WAS_IN_QUEUE)
  {
#ifdef ITER_DEBUG
    itercount++;
#endif
    if (u >= first_thru_node || u == start_node)
    {
      for (i = Va[u]; i < Va[u+1]; i++) /* all neighbours of u */
      {
        v = Ea[i];      /* v is adjacent to u */
        assert(v >= 0 && v < num_nodes);
        uvdist = Wa[i]; /* with this cost on edge */
        newcost = dist[u] + uvdist;
        if (newcost < dist[v])
        {
          dist[v] = newcost;
          prevlink[v] = i;
          if (queue_next[v] == WAS_IN_QUEUE)
          {
            queue_next[v] = queue_first;
            queue_first = v;
            if (queue_last == INVALID)
            {
              queue_last = v;
            }
          }
          else if (queue_next[v] == INVALID && v != queue_last)
          {
            if (queue_last != INVALID)
            {
              queue_next[queue_last] = v;
              queue_last = v;
            }
            else
            {
              queue_first = v;
              queue_last = v;
              queue_next[queue_last] = INVALID;
            }
          }
        }
      }
    }
    u = queue_first;
    if (u == INVALID || u == WAS_IN_QUEUE)
      break;
    assert(u >=0 && u < num_nodes);
    queue_first = queue_next[u];
    queue_next[u] = WAS_IN_QUEUE;
    if (queue_last == u)
      queue_last = INVALID;
  }
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &endtime);
  otime = 1000 * endtime.ru_utime.tv_sec + endtime.ru_utime.tv_usec/1000
          + 1000 * endtime.ru_stime.tv_sec + endtime.ru_stime.tv_usec/1000
          - 1000 * starttime.ru_utime.tv_sec - starttime.ru_utime.tv_usec/1000
          - 1000 * starttime.ru_stime.tv_sec - starttime.ru_stime.tv_usec/1000;
  printf("d'Esopo-Pape shortest path computatoin time: %d ms\n", otime);
#endif
#ifdef ITER_DEBUG
  printf("d'Esopo-Pape shortest path iterations: %d\n", itercount);
#endif

}


/*
 * get_shortest_path() - get shortest path from predecessor array
 *
 * Parameters:
 *      pred - predecessor array
 *      orig - source node
 *      dest - destination node
 *      path - (OUT) shortest path from source to dest NB in reverse
 *             must have enough space for path (max num nodes)
 *      pathlen (OUT) number of nodes in path (entries used in path array)
 *
 * Return value:
 *      None.
 */
void get_shortest_path(long pred[], long orig, long dest, long path[], 
                       long *pathlen)
{
  long v;

  v = dest;
  *pathlen = 0;
  if (pred[v] == INVALID)
    return; /* no path */
  while (v != orig)
  {
    assert(v != INVALID);
    path[(*pathlen)++] = v;
    v = pred[v];
  }
  path[(*pathlen)++] = orig;
}


/*
 * sssp_pape_lll() - single-source shortest path by d'Esopo-Pape algorithm
 *          
 *
 *  single-source shortest path by d'Esopo-Pape algorithm with LLL
 * (large label last) modification (see Klunder & Post 2006 and
 * citation [5] therein Bertsekas et al 1993 Networks 23:703-709).
 * In this modification of the algorithm, at each iteration,
 * when the node at the top f
 * the queue has a larger label (current distance) than the average
 * lable of nodes in the queue, this (top) node is not removed as current
 * node but instaed put at bottom of queue (and we keep doing this
 * until a node wih smaller than average label is found).
 *
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    start_node - source node
 *    first_thru_node - first node number that is allowed in a path
 *                      (earlier ones are actually 'zones' for origin/dest).
 *    prevlink  -  (OUT) predecessor link array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is link index of edge into node i
 *              in shortest path to node i
 *   dist (WORK) - vector of num_nodes doubles for distance from origin of each
 *   queue_next (WORK) - vector of num_nodes longs for queue of nodes
 */
void sssp_pape_lll(long Va[],long Ea[], double Wa[], 
                   long num_nodes, long num_edges,
                   long start_node, long first_thru_node, long prevlink[],
                   double dist[], long queue_next[])
{

  long u,v,next_u;
  long i;
  long queue_first, queue_last;
  double uvdist, newcost;
  double average_cost, total_cost; 
  long queue_size;
  double oldcost;
  int num_checked;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif
#ifdef SLOW_INTERNAL_DEBUG
  double tmp_total_cost;
  long queue_count;
  long w;
#endif
#ifdef ITER_DEBUG
  int itercount=0;
#endif

#ifdef USE_KAHAN_SUM
double kahan_c;
double kahan_y, kahan_t;
double kahan_test;
#define INIT_KAHAN_SUM(sum, init) {(sum) = (init); kahan_c = 0;  \
                                   kahan_test=(init);}
#define KAHAN_SUM(sum, summand) {kahan_y = (summand) - kahan_c;          \
                                 kahan_t = (sum) + kahan_y;              \
                                 kahan_c = (kahan_t - (sum)) - kahan_y;  \
                                 (sum) = kahan_t; \
                                 kahan_test += (summand);}
#else
#define INIT_KAHAN_SUM(sum,init) {(sum) = (init);}
//#define INIT_KAHAN_SUM(sum,init) {if ((init)>=FLOATINF)fprintf(stderr, "KAHAN_INIT >= FLOATINF\n"); (sum) = (init);}
#define KAHAN_SUM(sum, summand)  {(sum) += (summand);}
//#define KAHAN_SUM(sum, summand)  {if((summand)>=FLOATINF)fprintf(stderr, "KAHAN summand >= FLOATINF\n");(sum) += (summand);}
#endif /* USE_KAHAN_SUM */
                                 

  assert (!(start_node >= num_nodes));

#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prevlink[v] = INVALID;
    queue_next[v] = INVALID;
    if (v == start_node)
      dist[v] = 0;
    else
      dist[v] = FLOATINF;
  }
  queue_first = INVALID;
  queue_last = INVALID;
  queue_size = 0;
  INIT_KAHAN_SUM(total_cost, 0);

  u = start_node;
  while (u != INVALID && u != WAS_IN_QUEUE)
  {
#ifdef ITER_DEBUG
    itercount++;
#endif
    if (u >= first_thru_node || u == start_node)
    {
      for (i = Va[u]; i < Va[u+1]; i++) /* all neighbours of u */
      {
        v = Ea[i];      /* v is adjacent to u */
        assert(v >= 0 && v < num_nodes);
        uvdist = Wa[i]; /* with this cost on edge */
        newcost = dist[u] + uvdist;
        if (newcost < dist[v])
        {
          oldcost = dist[v];
          dist[v] = newcost;
          prevlink[v] = i;
          if (queue_next[v] == WAS_IN_QUEUE) /* if v was in queue */
          {                                
            queue_next[v] = queue_first;     /* add as first in queue */  
            queue_first = v;
            if (queue_last == INVALID)
            {
              queue_last = v;
            }
            ++queue_size;
            KAHAN_SUM(total_cost, dist[v]);
          }
          else if (queue_next[v]== INVALID && v != queue_last)/* wasn't in q*/
          { 
            /* if not in q and not in q before 
               add at end of queue */
            if (queue_last != INVALID) 
            {
              queue_next[queue_last] = v;
              queue_last = v;
              ++queue_size;
              KAHAN_SUM(total_cost, dist[v]);
            }
            else
            {
              queue_first = v;
              queue_last = v;
              queue_next[queue_last] = INVALID;
              queue_size = 1;
              INIT_KAHAN_SUM(total_cost, dist[v]);
            }
          }
          else   /* already in queue */
          {
//            total_cost -= oldcost - newcost; /* newcost < oldcost */
            //KAHAN_SUM(total_cost, newcost - oldcost); /* newcost < oldcost */
            KAHAN_SUM(total_cost, newcost);KAHAN_SUM(total_cost, -oldcost);
          }
        }
      }
    }

    /* LLL - check if first in queue is > average label in queue
       and if it is then put at end of queue, keep doing  so until
       we get head of queue wih <= average label as current node u
    */

#ifdef SLOW_INTERNAL_DEBUG    
    /* get average label in queue */
    tmp_total_cost = 0;
    queue_count = 0;
    for (w = queue_first; w != INVALID; w = queue_next[w])
    {
      if (dist[w] >= FLOATINF)
        fprintf(stderr, "dist[%d] = FLOATINF\n", w);
      tmp_total_cost += dist[w];
      ++queue_count;
    }
//    fprintf(stderr, "queue_size = %d queue_count = %d\n", queue_size, queue_count); // XXX
    assert(queue_count == queue_size);
    average_cost =  tmp_total_cost / queue_count;
//    fprintf(stderr, "total_cost = %f tmp_total_cost = %f\n", total_cost, tmp_total_cost); //XXX
 //   fprintf(stderr, "average_cost = %f total_cost/queue_size = %f\n", average_cost, total_cost/queue_size); //XXX
    //assert(fabs(total_cost - tmp_total_cost) < 1e-04); // large accumulation of fp rounding errors from adding differences
//    if(fabs(total_cost - tmp_total_cost) >= 1e-04)  // large accumulation of fp rounding errors from adding differences
 //   fprintf(stderr, "total_cost = %.15f tmp_total_cost = %.15f (%d in queue)\n", total_cost, tmp_total_cost, queue_size ); //XXX
//  fprintf(stderr, "avg cost (%d in queue) is %f\n", queue_size,average_cost);//XXX
#endif
    
    average_cost = total_cost / queue_size;

    num_checked = 0;
    u = queue_first;
    while (num_checked < queue_size)
    {
      assert(u >=0 && u < num_nodes);
      next_u = queue_next[u];
//XXX      fprintf(stderr, "cost of %d is %f\n", u, dist[u]);
      if ( queue_size > 1 &&
           dist[u] > average_cost &&
           num_checked < queue_size - 1)
      {
        /* label of u is larger than avg in queue, move it to end of queue */
        queue_first = queue_next[u];
        queue_next[u] = INVALID;
        if (queue_last != INVALID) 
        {
          queue_next[queue_last] = u;
          queue_last = u;
        }
        else
        {
          queue_first = u;
          queue_last = u;
        }
      }
      else
      {
        if ( ( queue_size > 1 && dist[u] > average_cost) &&
             !( num_checked < queue_size - 1 ) )
           fprintf(stderr, "numchecked=%d >= queue_size=%ld -1\n", num_checked, queue_size);//XXX

        /* label of u is not larger than average in queue, make it current */
        /* (also do this if only one left unchecked in queue */
        /*  or if we have checked all in queue and not found one larger than */
        /*  average: although this is "impossible", it actually can happen */
        /*  due to the accumulation of fp errors in summting the total_cost */
        /*  so the fast calcualted average_cost is not actually correct, as */
        /*  per the SLOW_INTERNAL_DEBUG code above. */
        queue_first = queue_next[u];
        queue_next[u] = WAS_IN_QUEUE;
        if (queue_last == u)
          queue_last = INVALID;
        --queue_size;
        KAHAN_SUM(total_cost, -dist[u]);
        break;
      }
      ++num_checked;
      u = next_u;
    }
  }
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &endtime);
  otime = 1000 * endtime.ru_utime.tv_sec + endtime.ru_utime.tv_usec/1000
          + 1000 * endtime.ru_stime.tv_sec + endtime.ru_stime.tv_usec/1000
          - 1000 * starttime.ru_utime.tv_sec - starttime.ru_utime.tv_usec/1000
          - 1000 * starttime.ru_stime.tv_sec - starttime.ru_stime.tv_usec/1000;
  printf("d'Esopo-Pape LLL shortest path computatoin time: %d ms\n", otime);
#endif
#ifdef ITER_DEBUG
  printf("d'Esopo-Pape LLL shortest path iterations: %d\n", itercount);
#endif


}


/*
 * sssp_slf() - single-source shortest path by small label first (SLF)
 *
 *   SLF algorithm: see Klunder & Post 2006 and citation [3] therein:
 *    Bertsekas 1993 "A Simple and Fast Label Correcting Algorithm for
 *     Shortest Paths" Networks 23:703-709
 *  Instead of enqueuing according to whether previously in queue or not,
 *  whenveer a node is enqueued, its label is compared with that at the top
 *  of the queue, and if its label is smaller than that at th etop,
 *  it is put at the top, otherwise the bottom of thq eueue.
 *  

 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    start_node - source node
 *    first_thru_node - first node number that is allowed in a path
 *                      (earlier ones are actually 'zones' for origin/dest).
 *    prevlink  -  (OUT) predecessor link array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is link index of edge into node i
 *              in shortest path to node i
 *   dist (WORK) - vector of num_nodes doubles for distance from origin of each
 *   queue_next (WORK) - vector of num_nodes longs for queue of nodes
 *    prevlink  -  (OUT) predecessor  array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is predecessor of node i
 *              in shortest path to node i
 */
void sssp_slf(long Va[],long Ea[], double Wa[],
                   long num_nodes, long num_edges,
                   long start_node, long first_thru_node, long prevlink[],
                   double dist[], long queue_next[])

{
  long u,v;
  long i;
  long queue_first, queue_last;
  float uvdist, newcost;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif


  assert (!(start_node >= num_nodes));

#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prevlink[v] = INVALID;
    queue_next[v] = INVALID;
    if (v == start_node)
      dist[v] = 0;
    else
      dist[v] = FLOATINF;
  }
  queue_first = INVALID;
  queue_last = INVALID;

  u = start_node;
  while (u != INVALID)
  {
    if (u >= first_thru_node || u == start_node)
    {
      for (i = Va[u]; i < Va[u+1]; i++) /* all neighbours of u */
      {
        v = Ea[i];      /* v is adjacent to u */
        assert(v >= 0 && v < num_nodes);
        uvdist = Wa[i]; /* with this cost on edge */
        newcost = dist[u] + uvdist;
        if (newcost < dist[v])
        {
          dist[v] = newcost;
          prevlink[v] = i;

          if (queue_next[v] == INVALID && v != queue_last)
          { 
            if (queue_first != INVALID && dist[v] < dist[queue_first])
            {     
              /* add as first in queue */  
              queue_next[v] = queue_first;
              queue_first = v;
              if (queue_last ==  INVALID)
                queue_last = queue_next[v];
            }
            else
            {
              /* add at end of queue */
              if (queue_last != INVALID && v != queue_last) 
              {
                assert(queue_last != v);
                queue_next[queue_last] = v;
                queue_last = v;
            assert(queue_next[queue_last] == INVALID);
              }
              else
              {
                if (queue_first == INVALID)
                {
                  queue_first = v;
                  queue_last = v;
                }
                else
                {
                  assert(queue_first != INVALID);
                  assert(queue_last == INVALID);
                  assert(queue_next[queue_first] == INVALID);
                  queue_next[queue_first] = v;
                  queue_last = v;
                }
              }
            }
            assert(queue_next[queue_last] == INVALID);
          }
        }
      }
    }
    /* get first node in queue as current node u */
    u = queue_first;
    if (u == INVALID)
      break;
    assert(u >=0 && u < num_nodes);
    queue_first = queue_next[u];
    queue_next[u] = INVALID;
    if (queue_last == u)
      queue_last = INVALID;
  }
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &endtime);
  otime = 1000 * endtime.ru_utime.tv_sec + endtime.ru_utime.tv_usec/1000
          + 1000 * endtime.ru_stime.tv_sec + endtime.ru_stime.tv_usec/1000
          - 1000 * starttime.ru_utime.tv_sec - starttime.ru_utime.tv_usec/1000
          - 1000 * starttime.ru_stime.tv_sec - starttime.ru_stime.tv_usec/1000;
  printf("SLF shortest path computatoin time: %d ms\n", otime);
#endif
}


/*
 * sssp_slf_lll() - single-source shortest path by small label first (SLF)
 *               with LLL modification
 *
 *   SLF algorithm: see Klunder & Post 2006 and citation [3] therein:
 *    Bertsekas 1993 "A Simple and Fast Label Correcting Algorithm for
 *     Shortest Paths" Networks 23:703-709
 *  Instead of enqueuing according to whether previously in queue or not,
 *  whenveer a node is enqueued, its label is compared with that at the top
 *  of the queue, and if its label is smaller than that at th etop,
 *  it is put at the top, otherwise the bottom of thq eueue.
 *  

 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    start_node - source node
 *    first_thru_node - first node number that is allowed in a path
 *                      (earlier ones are actually 'zones' for origin/dest).
 *    prevlink  -  (OUT) predecessor link array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is link index of edge into node i
 *              in shortest path to node i
 *   dist (WORK) - vector of num_nodes doubles for distance from origin of each
 *   queue_next (WORK) - vector of num_nodes longs for queue of nodes
 *    prevlink  -  (OUT) predecessor  array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is predecessor of node i
 *              in shortest path to node i
 */
void sssp_slf_lll(long Va[],long Ea[], double Wa[],
                   long num_nodes, long num_edges,
                   long start_node, long first_thru_node, long prevlink[],
                   double dist[], long queue_next[])

{

  long u,v;
  long next_u;
  long i;
  long queue_first, queue_last;
  double uvdist, newcost;
  double average_cost, total_cost;
  long queue_size;
  double oldcost;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif
#ifdef SLOW_INTERNAL_DEBUG
  double tmp_total_cost;
  long queue_count;
  long w;
#endif
  int num_checked;

  assert (!(start_node >= num_nodes));

#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prevlink[v] = INVALID;
    queue_next[v] = INVALID;
    if (v == start_node)
      dist[v] = 0;
    else
      dist[v] = FLOATINF;
  }
  queue_first = INVALID;
  queue_last = INVALID;
  queue_size = 0;
  total_cost = 0;

  u = start_node;
  while (u != INVALID)
  {
    if (u >= first_thru_node || u == start_node)
    {
      for (i = Va[u]; i < Va[u+1]; i++) /* all neighbours of u */
      {
        v = Ea[i];      /* v is adjacent to u */
        assert(v >= 0 && v < num_nodes);
        uvdist = Wa[i]; /* with this cost on edge */
        newcost = dist[u] + uvdist;
        if (newcost < dist[v])
        {
          oldcost = dist[v];
          dist[v] = newcost;
          prevlink[v] = i;

          if (queue_next[v] == INVALID && v != queue_last)
          { 
            if (queue_first != INVALID && dist[v] < dist[queue_first])
            {     
              /* add as first in queue */  
              queue_next[v] = queue_first;
              queue_first = v;
              if (queue_last ==  INVALID)
                queue_last = v;
              ++queue_size;
              total_cost += dist[v];
            }
            else
            {
              /* add at end of queue */
              if (queue_last != INVALID)
              {
                assert(queue_last != v);
                queue_next[queue_last] = v;
                queue_last = v;
                assert(queue_next[queue_last] == INVALID);
                ++queue_size;
                total_cost += dist[v];
              }
              else
              {
                queue_first = v;
                queue_last = v;
                queue_next[queue_last] = INVALID;
                queue_size = 1;
                total_cost = dist[v];
                /* if (queue_first == INVALID) */
                /* { */
                /*   queue_first = v; */
                /*   queue_last = v; */
                /*   queue_size = 1; */
                /*   total_cost = dist[v]; */
                /* } */
                /* else */
                /* { */
                /*   assert(queue_first != INVALID); */
                /*   assert(queue_last == INVALID); */
                /*   assert(queue_next[queue_first] == INVALID); */
                /*   queue_next[queue_first] = v; */
                /*   queue_last = v; */
                /*   queue_size = 2; */
                /*   total_cost = dist[v] + dist[queue_first]; */
                /* } */
              }
            }
            assert(queue_next[queue_last] == INVALID);
          }
          else /* already in queue */
          {
            total_cost += newcost - oldcost; /* newcost < oldcost */
          }
        }
      }
    }

    /* LLL - check if first in queue is > average label in queue
       and if it is then put at end of queue, keep doing  so until
       we get head of queue wih <= average label as current node u
    */

#ifdef SLOW_INTERNAL_DEBUG    
    /* get average label in queue */
    tmp_total_cost = 0;
    queue_count = 0;
    for (w = queue_first; w != INVALID; w = queue_next[w])
    {
      if (dist[w] >= FLOATINF)
        fprintf(stderr, "dist[%d] = FLOATINF\n", w);
      tmp_total_cost += dist[w];
      ++queue_count;
    }
//    fprintf(stderr, "queue_size = %d queue_count = %d\n", queue_size, queue_count); // XXX
    assert(queue_count == queue_size);
    average_cost =  tmp_total_cost / queue_count;
      fprintf(stderr, "total_cost = %f tmp_total_cost = %f\n", total_cost, tmp_total_cost); //XXX
   fprintf(stderr, "average_cost = %f total_cost/queue_size = %f\n", average_cost, total_cost/queue_size); //XXX
    assert(fabs(total_cost - tmp_total_cost) < 1e-04); // large accumulation of fp rounding errors from adding differences
//    if(fabs(total_cost - tmp_total_cost) >= 1e-04)  // large accumulation of fp rounding errors from adding differences
 //   fprintf(stderr, "total_cost = %.15f tmp_total_cost = %.15f (%d in queue)\n", total_cost, tmp_total_cost, queue_size ); //XXX
//  fprintf(stderr, "avg cost (%d in queue) is %f\n", queue_size,average_cost);//XXX
#endif
    
    average_cost = total_cost / queue_size;

    num_checked = 0;
    u = queue_first;
    while (num_checked < queue_size)
    {
      assert(u >=0 && u < num_nodes);
      next_u = queue_next[u];
//XXX      fprintf(stderr, "cost of %d is %f\n", u, dist[u]);
      if ( queue_size > 1 &&
           dist[u] > average_cost &&
           num_checked < queue_size - 1)
      {
        /* label of u is larger than avg in queue, move it to end of queue */
        queue_first = queue_next[u];
        queue_next[u] = INVALID;
        if (queue_last != INVALID) 
        {
          queue_next[queue_last] = u;
          queue_last = u;
        }
        else
        {
          queue_first = u;
          queue_last = u;
        }
      }
      else
      {
        if ( ( queue_size > 1 && dist[u] > average_cost) &&
             !( num_checked < queue_size - 1 ) )
           fprintf(stderr, "numchecked=%d >= queue_size=%ld -1\n", num_checked, queue_size);//XXX

        /* label of u is not larger than average in queue, make it current */
        /* (also do this if only one left unchecked in queue */
        /*  or if we have checked all in queue and not found one larger than */
        /*  average: although this is "impossible", it actually can happen */
        /*  due to the accumulation of fp errors in summting the total_cost */
        /*  so the fast calcualted average_cost is not actually correct, as */
        /*  per the SLOW_INTERNAL_DEBUG code above. */
        queue_first = queue_next[u];
        queue_next[u] = INVALID;
        if (queue_last == u)
          queue_last = INVALID;
        --queue_size;
        total_cost -= dist[u];
        break;
      }
      ++num_checked;
      u = next_u;
    }
  }
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &endtime);
  otime = 1000 * endtime.ru_utime.tv_sec + endtime.ru_utime.tv_usec/1000
          + 1000 * endtime.ru_stime.tv_sec + endtime.ru_stime.tv_usec/1000
          - 1000 * starttime.ru_utime.tv_sec - starttime.ru_utime.tv_usec/1000
          - 1000 * starttime.ru_stime.tv_sec - starttime.ru_stime.tv_usec/1000;
  printf("SLF_LLL shortest path computatoin time: %d ms\n", otime);
#endif


}








/*
 * sssp_prevnodefirst()_lll - single-source shortest path by previous node first
 *                            with LLL  modification.
 *
 *  Instead of enqueuing according to whether previously in queue or not,
 *  whenveer a node is enqueued, if it is the predecessor node for the current
 *  node from a previous run (input frmolast iteration) it is put on the
 * top of the queue otherwise the bottom.
 *  
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    start_node - source node
 *    first_thru_node - first node number that is allowed in a path
 *                      (earlier ones are actually 'zones' for origin/dest).
 *    old_prevnode - Each prevnode[i] is previous node in shortest path to i,
 *                   as set by this function in prevnode output vector.
 *                  On first call just set to INVALID in all.
 *    prevlink  -  (OUT) predecessor link array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is link index of edge into node i
 *              in shortest path to node i
 *    prevnode  -  (OUT) predecessor node array,
 *               must have space for num_nodes entries
 *              Each prevnode[i] is previous node in shortest path to i
 *   dist (WORK) - vector of num_nodes doubles for distance from origin of each
 *   queue_next (WORK) - vector of num_nodes longs for queue of nodes
 *    prevlink  -  (OUT) predecessor  array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is predecessor of node i
 *              in shortest path to node i
 */
void sssp_prevnodefirst_lll(long Va[],long Ea[], double Wa[],
                        long num_nodes, long num_edges,
                        long start_node, long first_thru_node, 
                        long old_prevnode[],
                        long prevlink[], long prevnode[],
                        double dist[], long queue_next[])
{
  long u,v;
  long next_u;
  long i;
  long queue_first, queue_last;
  double uvdist, newcost;
  double average_cost, total_cost;
  long queue_size;
  double oldcost;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif
#ifdef SLOW_INTERNAL_DEBUG
  double tmp_total_cost;
  long queue_count;
  long w;
#endif
  int num_checked;

  assert (!(start_node >= num_nodes));

#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prevlink[v] = INVALID;
    queue_next[v] = INVALID;
    if (v == start_node)
      dist[v] = 0;
    else
      dist[v] = FLOATINF;
  }
  queue_first = INVALID;
  queue_last = INVALID;
  queue_size = 0;
  total_cost = 0;

  u = start_node;
  while (u != INVALID)
  {
    if (u >= first_thru_node || u == start_node)
    {
      for (i = Va[u]; i < Va[u+1]; i++) /* all neighbours of u */
      {
        v = Ea[i];      /* v is adjacent to u */
        assert(v >= 0 && v < num_nodes);
        uvdist = Wa[i]; /* with this cost on edge */
        newcost = dist[u] + uvdist;
        if (newcost < dist[v])
        {
          oldcost = dist[v];
          dist[v] = newcost;
          prevlink[v] = i;
          prevnode[v] = u;

          if (queue_next[v] == INVALID && v != queue_last)
          { 
            if (queue_first != INVALID && u == old_prevnode[v])
            {     
              /* add as first in queue */  
              queue_next[v] = queue_first;
              queue_first = v;
              if (queue_last ==  INVALID)
                queue_last = v;
              ++queue_size;
              total_cost += dist[v];
            }
            else
            {
              /* add at end of queue */
              if (queue_last != INVALID)
              {
                assert(queue_last != v);
                queue_next[queue_last] = v;
                queue_last = v;
                assert(queue_next[queue_last] == INVALID);
                ++queue_size;
                total_cost += dist[v];
              }
              else
              {
                queue_first = v;
                queue_last = v;
                queue_next[queue_last] = INVALID;
                queue_size = 1;
                total_cost = dist[v];
                /* if (queue_first == INVALID) */
                /* { */
                /*   queue_first = v; */
                /*   queue_last = v; */
                /*   queue_size = 1; */
                /*   total_cost = dist[v]; */
                /* } */
                /* else */
                /* { */
                /*   assert(queue_first != INVALID); */
                /*   assert(queue_last == INVALID); */
                /*   assert(queue_next[queue_first] == INVALID); */
                /*   queue_next[queue_first] = v; */
                /*   queue_last = v; */
                /*   queue_size = 2; */
                /*   total_cost = dist[v] + dist[queue_first]; */
                /* } */
              }
            }
            assert(queue_next[queue_last] == INVALID);
          }
          else /* already in queue */
          {
            total_cost += newcost - oldcost; /* newcost < oldcost */
          }
        }
      }
    }

    /* LLL - check if first in queue is > average label in queue
       and if it is then put at end of queue, keep doing  so until
       we get head of queue wih <= average label as current node u
    */

#ifdef SLOW_INTERNAL_DEBUG    
    /* get average label in queue */
    tmp_total_cost = 0;
    queue_count = 0;
    for (w = queue_first; w != INVALID; w = queue_next[w])
    {
      if (dist[w] >= FLOATINF)
        fprintf(stderr, "dist[%d] = FLOATINF\n", w);
      tmp_total_cost += dist[w];
      ++queue_count;
    }
//    fprintf(stderr, "queue_size = %d queue_count = %d\n", queue_size, queue_count); // XXX
    assert(queue_count == queue_size);
    average_cost =  tmp_total_cost / queue_count;
      fprintf(stderr, "total_cost = %f tmp_total_cost = %f\n", total_cost, tmp_total_cost); //XXX
   fprintf(stderr, "average_cost = %f total_cost/queue_size = %f\n", average_cost, total_cost/queue_size); //XXX
    assert(fabs(total_cost - tmp_total_cost) < 1e-04); // large accumulation of fp rounding errors from adding differences
//    if(fabs(total_cost - tmp_total_cost) >= 1e-04)  // large accumulation of fp rounding errors from adding differences
 //   fprintf(stderr, "total_cost = %.15f tmp_total_cost = %.15f (%d in queue)\n", total_cost, tmp_total_cost, queue_size ); //XXX
//  fprintf(stderr, "avg cost (%d in queue) is %f\n", queue_size,average_cost);//XXX
#endif
    
    average_cost = total_cost / queue_size;

    num_checked = 0;
    u = queue_first;
    while (num_checked < queue_size)
    {
      assert(u >=0 && u < num_nodes);
      next_u = queue_next[u];
//XXX      fprintf(stderr, "cost of %d is %f\n", u, dist[u]);
      if ( queue_size > 1 &&
           dist[u] > average_cost &&
           num_checked < queue_size - 1)
      {
        /* label of u is larger than avg in queue, move it to end of queue */
        queue_first = queue_next[u];
        queue_next[u] = INVALID;
        if (queue_last != INVALID) 
        {
          queue_next[queue_last] = u;
          queue_last = u;
        }
        else
        {
          queue_first = u;
          queue_last = u;
        }
      }
      else
      {
        if ( ( queue_size > 1 && dist[u] > average_cost) &&
             !( num_checked < queue_size - 1 ) )
           fprintf(stderr, "numchecked=%d >= queue_size=%ld -1\n", num_checked, queue_size);//XXX

        /* label of u is not larger than average in queue, make it current */
        /* (also do this if only one left unchecked in queue */
        /*  or if we have checked all in queue and not found one larger than */
        /*  average: although this is "impossible", it actually can happen */
        /*  due to the accumulation of fp errors in summting the total_cost */
        /*  so the fast calcualted average_cost is not actually correct, as */
        /*  per the SLOW_INTERNAL_DEBUG code above. */
        queue_first = queue_next[u];
        queue_next[u] = INVALID;
        if (queue_last == u)
          queue_last = INVALID;
        --queue_size;
        total_cost -= dist[u];
        break;
      }
      ++num_checked;
      u = next_u;
    }
  }


#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &endtime);
  otime = 1000 * endtime.ru_utime.tv_sec + endtime.ru_utime.tv_usec/1000
          + 1000 * endtime.ru_stime.tv_sec + endtime.ru_stime.tv_usec/1000
          - 1000 * starttime.ru_utime.tv_sec - starttime.ru_utime.tv_usec/1000
          - 1000 * starttime.ru_stime.tv_sec - starttime.ru_stime.tv_usec/1000;
  printf("sssp_prevnodefirst shortest path computatoin time: %d ms\n", otime);
#endif
}

