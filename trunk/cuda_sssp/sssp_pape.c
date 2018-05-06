/*****************************************************************************
 * 
 * File:    sssp_pape.c
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: sssp_pape.c 435 2011-06-29 02:07:20Z astivala $
 *
 * single-source shortest paths using 
 * the d'Esopo-Pape algorithm [Pape 1974 "Implementation and Efficiency
 * of Moore-Algorithms for the shortest route problem" Math. Prorgam. 7:212-222]
 * and variations of it.
 * This algorithm is well suited to shorest paths in road networks
 * (see Klunder & Post 2006 "The Shortest Path Problem on Large-Scale
 * Real-Road Networks" Networks 48(4):182-194).
 *
 * The d'Esopo-Pape is termed a "label correcting" rather than "label setting"
 * (Dijkstra type) algorithm. A node is always removed from head of queue,
 * and is placed at the tail of the queue if it has never been in the queue
 * before, otehrwised placed at the head.
 *
 * TODO bidrectional versions (see Klunder & Post 2006).
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

/****************************************************************************
 *
 * constants and type definitions
 *
 ****************************************************************************/

#define INVALID      -1
#define WAS_IN_QUEUE -2


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
void dump_packed_arrays(const int Va[], const int Ea[], const float Wa[],
                        int num_nodes, int num_edges)
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
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    start_node - source node
 *    first_thru_node - first node number that is allowed in a path
 *                      (earlier ones are actually 'zones' for origin/dest).
 *    dist - (OUT) distance array, dist from source to each node
 *    prev  -  (OUT) predecessor  array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is predecessor of node i
 *              in shortest path to node i
 */

void sssp_pape(int Va[],int Ea[], float Wa[], int num_nodes, int num_edges,
               int start_node, int first_thru_node, float dist[], int prev[])
{

  int u,v;
  int i;
  int *queue_next;
  int queue_first, queue_last;
  float uvdist, newcost;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif


  assert (!(start_node >= num_nodes));


  if (!(queue_next = (int *)malloc(num_nodes *sizeof(int))))
  {
    fprintf(stderr, "malloc queue failed\n");
    return;
  }
  
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prev[v] = INVALID;
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
          prev[v] = u;
          if (queue_next[v] == WAS_IN_QUEUE) /* if v was in queue */
          {                                
            queue_next[v] = queue_first;     /* add as first in queue */  
            queue_first = v;
            if (queue_last == INVALID)
            {
              queue_last = v;
            }
          }
          else if (queue_next[v] == INVALID && v != queue_last)
          {
            if (queue_last != INVALID) /* if not in q and not in q before */
            {
              queue_next[queue_last] = v; /* add at end of queue */
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
    /* get first node in queue as current node u */
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

  free(queue_next);

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
void get_shortest_path(int pred[], int orig, int dest, int path[], 
                       int *pathlen)
{
  int v;

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
 *    dist - (OUT) distance array, dist from source to each node
 *    prev  -  (OUT) predecessor  array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is predecessor of node i
 *              in shortest path to node i
 */

void sssp_slf(int Va[],int Ea[], float Wa[], int num_nodes, int num_edges,
               int start_node, int first_thru_node, float dist[], int prev[])
{

  int u,v;
  int i;
  int *queue_next;
  int queue_first, queue_last;
  float uvdist, newcost;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif


  assert (!(start_node >= num_nodes));


  if (!(queue_next = (int *)malloc(num_nodes *sizeof(int))))
  {
    fprintf(stderr, "malloc queue failed\n");
    return;
  }
  
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prev[v] = INVALID;
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
          prev[v] = u;

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

  free(queue_next);

}


/*
 * sssp_pape_lll() - single-source shortest path by d'Esopo-Pape algorithm
 *          
 *
 *  single-source shortest path by d'Esopo-Pape algorithm with LLL
 * (large label last) modification (see Klunder & Post 2006 and
 * citation [5] therein Bertsekas et al 1996 J. Optimization Theory and
 * Applications 88(2):297-320.
 * In this modification of the algorithm, at each iteration,
 * when the node at the top of
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
 *    dist - (OUT) distance array, dist from source to each node
 *    prev  -  (OUT) predecessor  array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is predecessor of node i
 *              in shortest path to node i
 */

void sssp_pape_lll(int Va[],int Ea[], float Wa[],
                   int num_nodes, int num_edges,
                   int start_node, int first_thru_node, float dist[], int prev[])
{

  int u,v,w,next_u;
  int i;
  int *queue_next;
  int queue_first, queue_last;
  float uvdist, newcost;
  double average_cost, total_cost; // NB double to help avoid rounding problems
  float tmp_total_cost;
  int queue_size;
  int queue_count;
  float oldcost;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif


  assert (!(start_node >= num_nodes));


  if (!(queue_next = (int *)malloc(num_nodes *sizeof(int))))
  {
    fprintf(stderr, "malloc queue failed\n");
    return;
  }
  
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prev[v] = INVALID;
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
  while (u != INVALID && u != WAS_IN_QUEUE)
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
          prev[v] = u;
          if (queue_next[v] == WAS_IN_QUEUE) /* if v was in queue */
          {                                
            queue_next[v] = queue_first;     /* add as first in queue */  
            queue_first = v;
            if (queue_last == INVALID)
            {
              queue_last = v;
            }
            ++queue_size;
            total_cost += dist[v];
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
              total_cost += dist[v];
            }
            else
            {
              queue_first = v;
              queue_last = v;
              queue_next[queue_last] = INVALID;
              queue_size = 1;
              total_cost = dist[v];
            }
          }
          else   /* already in queue */
          {
//            total_cost -= oldcost - newcost; /* newcost < oldcost */
            total_cost += (double)newcost - (double)oldcost; /* newcost < oldcost */
          }
        }
      }
    }

    /* LLL - check if first in queue is > average label in queue
       and if it is then put at end of queue, keep doing  so until
       we get head of queue wih <= average label as current node u
    */

#undef SLOW_INTERNAL_DEBUG
#ifdef SLOW_INTERNAL_DEBUG    
    /* get average label in queue */
    tmp_total_cost = 0;
    queue_count = 0;
    for (w = queue_first; w != INVALID; w = queue_next[w])
    {
      tmp_total_cost += dist[w];
      ++queue_count;
    }
//   fprintf(stderr, "queue_size = %d queue_count = %d\n", queue_size, queue_count); // XXX
    assert(queue_count == queue_size);
    average_cost =  tmp_total_cost / queue_count;
    fprintf(stderr, "total_cost = %f tmp_total_cost = %f\n", total_cost, tmp_total_cost); //XXX
    fprintf(stderr, "average_cost = %f total_cost/queue_size = %f\n", average_cost, total_cost/queue_size); //XXX
    assert(fabs(total_cost - tmp_total_cost) < 0.1); // large accumulation of fp rounding errors from adding differences
//    fprintf(stderr, "avg cost (%d in queue) is %f\n", queue_size,average_cost);//XXX
#endif
    
    average_cost = total_cost / queue_size;

    u = queue_first;
    while (u != INVALID && u != WAS_IN_QUEUE)
    {
      assert(u >=0 && u < num_nodes);
      next_u = queue_next[u];
//XXX      fprintf(stderr, "cost of %d is %f\n", u, dist[u]);
      /* FIXME fix this to handle "impossible" situations due to fp rounding
         as per sssp_pape_lll() in tap/sssp_pape.c */

      if (queue_size > 1 && dist[u] > average_cost)
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
        /* label of u is not larger than average in queue, make it current */
        queue_first = queue_next[u];
        queue_next[u] = WAS_IN_QUEUE;
        if (queue_last == u)
          queue_last = INVALID;
        --queue_size;
        total_cost -= dist[u];
        break;
      }
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

  free(queue_next);

}


/*
 * sssp_bellmanford() - single-source shortest path by Bellman-Ford
 *  
 * Bellman-Ford algorithm: node is always added at bottom of queue.
 *
 * Parameters:
 *    Va, Ea, Wa - graph in packed adjacency list represention
 *    num_nodes - number of nodes (elemnts in  Va)
 *    num_edges - number of edges (elements in Ea, Wa)
 *    start_node - source node
 *    first_thru_node - first node number that is allowed in a path
 *                      (earlier ones are actually 'zones' for origin/dest).
 *    dist - (OUT) distance array, dist from source to each node
 *    prev  -  (OUT) predecessor  array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is predecessor of node i
 *              in shortest path to node i
 */

void sssp_bellmanford(int Va[],int Ea[], float Wa[], int num_nodes, int num_edges,
               int start_node, int first_thru_node, float dist[], int prev[])
{

  int u,v;
  int i;
  int *queue_next;
  int queue_first, queue_last;
  float uvdist, newcost;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif


  assert (!(start_node >= num_nodes));


  if (!(queue_next = (int *)malloc(num_nodes *sizeof(int))))
  {
    fprintf(stderr, "malloc queue failed\n");
    return;
  }
  
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prev[v] = INVALID;
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
          prev[v] = u;
          /* always at end is Bellman-Ford */
          if (queue_next[v] == INVALID)
          { 
            /* add at end of queue */
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
  printf("Bellman-Ford shortest path computatoin time: %d ms\n", otime);
#endif

  free(queue_next);

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
 *    dist - (OUT) distance array, dist from source to each node
 *    prev  -  (OUT) predecessor  array,
 *               must have space for num_nodes entries
 *              Each prevlink[i] is predecessor of node i
 *              in shortest path to node i
 */

void sssp_slf_lll(int Va[],int Ea[], float Wa[], int num_nodes, int num_edges,
               int start_node, int first_thru_node, float dist[], int prev[])
{

  int u,v;
  int w,next_u;
  int i;
  int *queue_next;
  int queue_first, queue_last;
  float uvdist, newcost;
  double average_cost, total_cost; // NB double to help avoid rounding problems
  float tmp_total_cost;
  int queue_size;
  int queue_count;
  float oldcost;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif


  assert (!(start_node >= num_nodes));


  if (!(queue_next = (int *)malloc(num_nodes *sizeof(int))))
  {
    fprintf(stderr, "malloc queue failed\n");
    return;
  }
  
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  for (v = 0; v < num_nodes; v++)
  {
    prev[v] = INVALID;
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
          prev[v] = u;

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
            queue_size++;
            total_cost += dist[v];
          }
          else /* already in queue */
          {
            total_cost += (double)newcost - (double)oldcost; /* newcost < oldcost */
          }
        }
      }
    }

    /* LLL - check if first in queue is > average label in queue
       and if it is then put at end of queue, keep doing  so until
       we get head of queue wih <= average label as current node u
    */

#undef SLOW_INTERNAL_DEBUG_SLF_LLL
#ifdef SLOW_INTERNAL_DEBUG_SLF_LLL
    /* get average label in queue */
    tmp_total_cost = 0;
    queue_count = 0;
    for (w = queue_first; w != INVALID; w = queue_next[w])
    {
      tmp_total_cost += dist[w];
      ++queue_count;
    }
//   fprintf(stderr, "queue_size = %d queue_count = %d\n", queue_size, queue_count); // XXX
    assert(queue_count == queue_size);
    average_cost =  tmp_total_cost / queue_count;
    fprintf(stderr, "total_cost = %f tmp_total_cost = %f\n", total_cost, tmp_total_cost); //XXX
    fprintf(stderr, "average_cost = %f total_cost/queue_size = %f\n", average_cost, total_cost/queue_size); //XXX
    assert(fabs(total_cost - tmp_total_cost) < 0.1); // large accumulation of fp rounding errors from adding differences
//    fprintf(stderr, "avg cost (%d in queue) is %f\n", queue_size,average_cost);//XXX
#endif
    
    average_cost = total_cost / queue_size;

    u = queue_first;
    while (u != INVALID)
    {
      assert(u >=0 && u < num_nodes);
      next_u = queue_next[u];
//XXX      fprintf(stderr, "cost of %d is %f\n", u, dist[u]);
      /* FIXME fix this to handle "impossible" situations due to fp rounding
         as per sssp_pape_lll() in tap/sssp_pape.c */
      if (queue_size > 1 && dist[u] > average_cost)
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
        /* label of u is not larger than average in queue, make it current */
        queue_first = queue_next[u];
        queue_next[u] = INVALID;
        if (queue_last == u)
          queue_last = INVALID;
        --queue_size;
        total_cost -= dist[u];
        break;
      }
      u = next_u;
    }
  }

#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &endtime);
  otime = 1000 * endtime.ru_utime.tv_sec + endtime.ru_utime.tv_usec/1000
          + 1000 * endtime.ru_stime.tv_sec + endtime.ru_stime.tv_usec/1000
          - 1000 * starttime.ru_utime.tv_sec - starttime.ru_utime.tv_usec/1000
          - 1000 * starttime.ru_stime.tv_sec - starttime.ru_stime.tv_usec/1000;
  printf("SLF shortest path computatoin time: %d ms\n", otime);
#endif

  free(queue_next);

}

