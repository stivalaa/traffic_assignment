/*****************************************************************************
 * 
 * File:    tap_frankwolfe_pthread.c
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_frankwolfe_pthread.c 823 2011-10-31 04:35:15Z astivala $
 *
 * Traffic assignment by Frank-Wolfe algorithm.
 * This version runs the shortest path computations (d'Esopo-Pape 
 * algorithm) in parallel using POSIX threads (one thread per
 * source node, up to max number of threads speicfieid to use).
 *
 * There are many descriptions of this algorithm, a simple one for
 * example is in Sheffi, Y. 1985 "Urban Transportation Networks:
 * Equilibribum Analysis with Mathematical Programming Methods"
 * Prentice-Hall (out of print, available online from
 * http://web.mit.edu/sheffi/www/urbanTransportation.html
 *
 *   Usage: tap_frankwolfe_pthread [-n threads] [-w flowsfilename] 
 *                                 [-i iter] [-r relgap] 
 *                                 netfilename  demandfilename
 *
 *     netfilename    is name of the net file defining node and links
 *     demandfilename is name of Origin-Destination demand file
 *     -n threads       : run with that many worker threads
 *     -w flowsfilename : warmstart from flows stored in flowsfilename
 *                        which is the output link flows file format
 *     -i iterations    : terminate after this many iterations
 *     -r relgap        : terminate after relative gap is <= relgap.
 *
 *     Output is link flows on stdout.
 *     
 *   Example usage:
 *   
 *   tap_frankwolfe_pthread SiouxFalls_net.txt SiouxFalls_trips.txt > SioxFalls_flow.txt
 *   
 *
 * These signals are handled:
 *
 *     SIGUSR1 (signal 10) - write status (current values, iteration etc.)
 *                           to stderr and continue. (UNLESS MPICH is used
 *                           as it needs SIGUSR1 itself)
 *
 *
 *     SIGUSR2 (signal 12) - write status (current values, iteration etc.)
 *                           to stderr and also write flow output file.
 *                           The -w option can be used to continue later
 *                           as a kind of checkpoint/restart.
 *
 *     SIGTERM (signal 15),
 *     SIGINT  (signal 2),
 *     SIGHUP  (signal 1)  
 *     SIGXCPU
 *     SIGPWR              - terminate cleanly at end of current iteration,
 *                           writing flow output file etc. just as if
 *                           the iteration or relative gap limit was reached.
 *                           The -w option can be used to continue later.
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <unistd.h>
#include <signal.h>
#include "parsetapfiles.h"
#include "sssp_pape.h"
#include "tap_functions.h"
#include "atomicdefs.h"
#include "utils.h"
#include "tap_frankwolfe_pthread.h"

#define TIME_SHORTEST_PATHS

/* TODO: create threads only once at beginning and signal them to 
 *       do new work rather than create/join each iteration
 *
 * TODO: paralleize various other loops (link volume idfference, update etc.)
 *       as marked in comments at appropriate points.
 */

/*
   In this program, all data for links will be in (1d) arrays such
   that the entries are sorted by link origin (ascending) and then by
   link destinatino (ascending). This way we can have net.links[k]
   being the data about link k and link_volumes[k] and link_costs[k]
   the current volume and cost for link k, the latter can be used
   directly as the weight Wa[k] in the packed adjacency list
   representation of the network (used for efficient shortest paths
   computation) since that is sorted this way.

   Note that node (and zone) numbers start at 1
   (since nodes are numbered in input and output files this way) but
   links start normal C style at the 0th entry in the array as they
   are not numbered externally.  
*/
   
/*****************************************************************************
 *
 * global data
 *
 ****************************************************************************/

unsigned int num_cores;   /* number of cores found on system */
unsigned int num_threads; /* maximum number of worker threads allowed */


   
/*****************************************************************************
 *
 * pthread data
 *
 ****************************************************************************/

/* thread data struct for the sssp_volume_update_thread() function */
typedef struct thread_data_s
{
  int thread_id; /* id of this thread to index into per-thread arrays */
                 /* 0..num_threads-1 NB NOT the pthread id */
  long *Va;
  long *Ea;
  double *Wa;
  long num_zones;
  long num_nodes;
  long num_edges;
  long first_start_node;   /* different for each thread */
  long last_start_node;    /* ierate first..last start nodes (inclusive) */
  long first_thru_node;
  long *prevlink;   /* different vector for each thread */
  double *link_volumes;
  demand_data_t **demands;
  link_data_t *links;
  double *dist; /* sssp_pape() workspace, different vector for each thread */
  long *queue_next; /* sssp_pape() workspace, differet vector for each thread */
} thread_data_t;

/* structures passed to each thread  as parameter */
static thread_data_t thread_data[MAX_NUM_THREADS];

/* pthreads thread handles */
static pthread_t threads[MAX_NUM_THREADS];


/* We will create a thread specifically to handle signals (just waits, does 
   nothing else at all), it will set global signal_status and main loop
   in tap_frankwolfe() will check it and perform appropriate action */

#define SIGNAL_STATUS_NONE       0 /* nothing; just keep going */
#define SIGNAL_STATUS_WRITEMSG   1 /* write status to stderr and continue */
#define SIGNAL_STATUS_CHECKPOINT 2 /* write status and flows file and conintue */
#define SIGNAL_STATUS_TERMINATE  3 /* write output and terminate cleanly */

static volatile sig_atomic_t signal_status = SIGNAL_STATUS_NONE;


/*****************************************************************************
 *
 * Local functions
 *
 ****************************************************************************/

/* 
 * atomic_add_double() - atomic add on double precision floating point
 *
 * Parameters:
 *     address - address of double to do min on
 *     val     - value to add to value at that address
 *
 * Return value:
 *     Old value
 *
 *  Read double at address as old and store old+val at address
 *  atomically.
 *
 */
static double atomic_add_double(double *address, double val)
{
  union double_longlong_union {
      double             dbl;
      unsigned long long ull;
  };
  union double_longlong_union old;
  union double_longlong_union assumed;
  union double_longlong_union newvalue;
  old.dbl = *address;
  do
  {
    assumed = old;
    newvalue.dbl = old.dbl + val;
    old.ull = CAS64((unsigned long long*)address, assumed.ull, newvalue.ull);
  }  
  while (assumed.ull != old.ull);
  return old.dbl;
}

/*
 * pthreads function to handle signals - does nothing else but wait for
 * wignal and set global  signal_status appropriately
 *
 * Parameters:
 *     threadarg - sigset_t for signals to handle
 *
 * Return value:
 *     None
 *
 * Uses global data:
 *    signal_status - set according to signal received
 */
static void *signal_handler_thread(void *threadarg)
{
  sigset_t *set = (sigset_t *)threadarg;
  int rc, sig;
  while (1)
  {
    if ((rc = sigwait(set, &sig) != 0))
    {
      perror("ERROR: sigwait() in signal handler thread failed, "
             "signals not handled");
      return NULL;
    }
    fprintf(stderr, "signal %d received...\n", sig);
    switch(sig)
    {
      case SIGUSR1:
        signal_status = SIGNAL_STATUS_WRITEMSG;
        break;

      case SIGUSR2:
        signal_status = SIGNAL_STATUS_CHECKPOINT;
        break;

      case SIGHUP:
      case SIGINT:
      case SIGTERM:
      case SIGXCPU:
      case SIGPWR:
        signal_status = SIGNAL_STATUS_TERMINATE;
        break;

      default:
        fprintf(stderr, "ERROR: unhandled signal %d\n", sig);
        break;
    }
    if (signal_status == SIGNAL_STATUS_TERMINATE)
      break;
  }
  return NULL;
}


/*
 * pthreads function to run sssp_pape() and update volume according
 * to resulting shortest path for a range of origin nodes
 * 
 */
static void *sssp_volume_update_thread(void *threadarg)
{
  thread_data_t *mydata = (thread_data_t *)threadarg;
  long orig;
  demand_data_t **demands = mydata->demands;
  long v,k,dest,i;
  double routeflow;
  int fail;

  for (orig = mydata->first_start_node; 
       orig <= mydata->last_start_node && orig < mydata->num_zones + 1;
       orig++)
  {
    /* get shortest paths to all nodes from this orgin */
    sssp_pape_lll(mydata->Va, mydata->Ea, mydata->Wa, mydata->num_nodes,
                  mydata->num_edges, orig, mydata->first_thru_node, 
                  mydata->prevlink, mydata->dist, mydata->queue_next);
    
    /* now process the results for this origin (start node) only */
    for (i = 0; (dest = demands[orig][i].dest) != 0; i++)
    {
      routeflow = demands[orig][i].demand;
      if (orig == dest || routeflow <= 0)/*should never be -ve; stops warning*/
        continue;
      v = dest;
      fail = 0;
      while (v != orig)
      {
        k = mydata->prevlink[v];
        if (k < 0)
        {
          /* this shouldn't happen, it means there is a nonzero demand
             between two nodes but no path between them  - 
             probably due to bad input data */
          fprintf(stderr, "ERROR: no path from %ld to %ld\n", orig, dest);
          fail = 1;
          exit(1);
        }
        assert(k >= 0);
        assert(mydata->links[k].term_node == v);
        /* link_volumes[k] += routeflow must be atomic as multiple of these
           threads may be updating the same link_volumes[k] simultaneously */
        atomic_add_double(&mydata->link_volumes[k], routeflow);
        v = mydata->links[k].init_node;
      }
    }
  }
  return NULL;
}




/*
 *  fw_assign() - flow assignment in Frank-Wolfe
 *
 *   Assign flow from O-D demands to links according to current shortest
 *   paths (lowest costs) on netgraph
 *   This version runs num_threads (module static var) threads 
 *   in parallel doing sssp_pape() for different origins
 *
 * Parameters:
 *    net - net structure parsed from net file
 *    demands  - demands array parsed from trips file
 *                for each origin demands[origin] is an array of demands structs
 *                terminated by one with 0  dest (origins and dests start at 1)
 *    Va, Ea, Wa - network in packed adjancey list format (Wa is current costs)
 *    link_volumes - (OUT) - volume on each link
 *   dist (WORK) - space for, for each of num_threads threads,
 *                  vector of num_nodes doubles for distance from origin of each
 *                 for sssp_pape()
 *   queue_next (WORK) - space for, for each of num_thread threads,
 *                       vector of num_nodes longs for queue of nodes
 *                 for sssp_pape()
 *   predlink (WORK) - space for, for each of num_threads threads,
 *                     vector of num_nodes longs for predecessor of each node
 *
 * Uses global data:
 *    num_threads - number of threads to use
 *
 *  Return value;
 *    0 if OK, nonzero on error.
 */
static int fw_assign(net_data_t *net, demand_data_t *demands[],
                      long Va[], long Ea[], double Wa[], 
                      double link_volumes[],
                      double dist[], long queue_next[], 
                      long predlink[])
{
  int  rc;
  unsigned int last_threadnum=0;
  unsigned int threadnum;
  long k;
#ifdef TIME_SHORTEST_PATHS
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int etime  = 0;
#endif
  long num_origins_per_thread;

  /* allocate space for a predlink[] vector for each thread */
  if (!(predlink = (long *)malloc(num_threads * 
                           (net->num_nodes+1)*sizeof(long))))
  {
    fprintf(stderr, "malloc pred failed\n");
    return(1);
  }

  for (k = 0; k < net->num_links; k++)
    link_volumes[k] = 0;

#ifdef TIME_SHORTEST_PATHS
  gettimeofday(&start_timeval, NULL);
#endif
  
  /* each thread does ceil(number_of_start_nodes, num_threads) origins */
  num_origins_per_thread = iDivUp(net->num_zones, num_threads);

  /* launch threads eaching doing num_origins_per_thread sssp + volume update */
  for (threadnum = 0; threadnum < num_threads; threadnum++)
  {
    thread_data[threadnum].thread_id = threadnum;
    thread_data[threadnum].Va = Va;
    thread_data[threadnum].Ea = Ea;
    thread_data[threadnum].Wa = Wa;
    thread_data[threadnum].num_zones = net->num_zones;
    thread_data[threadnum].num_nodes = net->num_nodes + 1;
    thread_data[threadnum].num_edges = net->num_links;
    thread_data[threadnum].first_start_node = 1 + 
                                             threadnum * num_origins_per_thread;
    thread_data[threadnum].last_start_node =  
      thread_data[threadnum].first_start_node + num_origins_per_thread - 1;
    thread_data[threadnum].first_thru_node = net->first_thru_node;
    thread_data[threadnum].prevlink = &predlink[threadnum*(net->num_nodes+1)];
    thread_data[threadnum].link_volumes = link_volumes;
    thread_data[threadnum].demands = demands;
    thread_data[threadnum].links = net->links;
    thread_data[threadnum].dist = &dist[threadnum*(net->num_nodes+1)];
    thread_data[threadnum].queue_next= &queue_next[threadnum*(net->num_nodes+1)];
    
//      fprintf(stderr, "starting thread %d\n", threadnum); // XXX
    if ((rc = pthread_create(&threads[threadnum], NULL,
                             sssp_volume_update_thread,
                             (void *)&thread_data[threadnum])))
    {
      fprintf(stderr, "pthread_create() failed (%d)\n", rc);
      return(1);
    }
    last_threadnum = threadnum;
  }
  
  /* wait for all the threads to finish */
//     fprintf(stderr, "last_threadnum = %d\n", last_threadnum); //XXX
  for (threadnum = 0; threadnum < num_threads; threadnum++)
  {
//      fprintf(stderr, "joining thread %d\n", threadnum); // XXX
    if ((rc = pthread_join(threads[threadnum], NULL)))
    {
      fprintf(stderr, "pthread_join () failed (%d)\n", rc);
      return(1);
    }
  }
  
#ifdef TIME_SHORTEST_PATHS
  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  etime = 1000*elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
  fprintf(stderr, "threaded shortest path + volume update total time %d ms\n", etime);
#endif

  free(predlink);
  return 0;
}


/* 
 * linkflow_entry_compar() - qsort comparison function for link_flow entries
 * 
 * Compares by 'from' node number first  then by 'to' node number if equal
 *
 */
static int linkflow_entry_compar(const void *ent1, const void *ent2)
{
  const link_flow_t *e1 = (const link_flow_t *)ent1;
  const link_flow_t *e2 = (const link_flow_t *)ent2;
  
  if (e1->init_node < e2->init_node)
    return -1;
  else if(e1->init_node > e2->init_node)
    return 1;
  else
    return ( e1->term_node < e2->term_node ? -1 : 
             (e1->term_node > e2->term_node ? 1 : 0) );
}


/*****************************************************************************
 *
 * Externally visible functions
 *
 ****************************************************************************/



/*
 * tap_frankwolfe() - solve traffic assignment problem by Frank-Wolfe method 
 * 
 * Parameters:
 *    net -  network (link) data parsed from net file
 *    demands - Origin-Dest demand data
 *    warm_start_mode - nonzero to read flows from previous run and start there
 *    flows_input_fp - only if warmstart_mode!=0, open (read) fp for flows input
 *    target_iterations - max number of iterations to run or 0
 *    target_relgap - stop when relative_gap <= relgap if target_relgap > 0
 *    flow_fp - open (write) filepointer to write flows output file
 *    total_cost_at_end (out) - total cost (VHT) 
 *    objvalue_at_end (out) - objective value at end
 *    iterations_at_end (out) - number of iterations run
 *    relgap_at_end (out) - relative gap reached
 *    time_each_iteration - if nonzero, output stats at each iteration to stderr
 *
 * Uses global data:
 *    num_threads - number of threads to use
 *
 * Return value:
 *    total cost at end
 */
int tap_frankwolfe(net_data_t *net,
                   demand_data_t **demands,
                    int warm_start_mode,
                    FILE *flows_input_fp,
                    int target_iterations,
                    double target_relgap,
                      FILE *flow_fp,
		   double *total_cost_at_end, double *objvalue_at_end,
		   int *iterations_at_end, double *relgap_at_end,
           int time_each_iteration)
{

  long *Va,*Ea;  /* vertex and edge arrays for packed adjlist format of graph */
  int i;
  double *link_volumes, *link_costs, *link_volumes2;
  int iter;
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int etime;
  double min_nonzero_capacity;
  double lambda;
  double total_demand = 0;
  double gap, relative_gap, lowerbound, best_lowerbound, objective_value;
  double link_cost_total, average_excess_cost;
  long flows_num_links = 0;
  link_flow_t *link_flows = NULL;
  double *dist;
  long *queue_next, *predlink;
  double total_cost;
  time_t now;
  int finished = 0;
  sigset_t set;
  int rc;
  pthread_t sigthread;

  /* block SIGUSR1, SIGUSR2, SIGHUP,SIGINT,SIGTERM, SIGXCPU, SIGPWR
     and set up signal handler thread to 
     handle them. Other threads crated will inherit the signal mask and
     also ignore them -only signal_handler_thread will deal with them */
  sigemptyset(&set);
#ifndef MPICH2
  /* MPICH needs SIGUSR1 for itself */
  sigaddset(&set, SIGUSR1);
#endif 
  sigaddset(&set, SIGUSR2);
  sigaddset(&set, SIGHUP);
  sigaddset(&set, SIGINT);
  sigaddset(&set, SIGTERM);
  sigaddset(&set, SIGXCPU);
  sigaddset(&set, SIGPWR);
  if ((rc = pthread_sigmask(SIG_BLOCK, &set, NULL)) != 0)
  {
    perror("ERROR: pthread_sigmask() failed");
    fprintf(stderr, "WARNING: signals will not be handled");
  }
  if ((rc = pthread_create(&sigthread, NULL, signal_handler_thread, &set)) != 0)
  {
    fprintf(stderr, "ERROR: create signal handler thread failed (%d)\n", rc);
    fprintf(stderr, "WARNING: signals will not be handled");
  }
                             

  if (warm_start_mode)
  {
    /* we ahve the 'warm start' option to parse link volumes (and costs)
       from a flows file (output of this program also) */
    assert(flows_input_fp != NULL);
    if (parse_flows_file(flows_input_fp, &link_flows, &flows_num_links) != 0)
    {
      fprintf(stderr, "error parsing flows input file\n");
      return(1);
    }
    fclose(flows_input_fp);

    /* sort by 'from' node ascending and within that by 'to' node ascending */
    qsort(link_flows,flows_num_links,sizeof(link_flow_t),linkflow_entry_compar);

    if (flows_num_links < net->num_links)
    {
      /* int num_new_links = net->num_links - flows_num_links; */
      int istart = 0, inserted_count = 0  ;
      /*links have been added in the net file relative to the flows input file*/
      if (!(link_flows = (link_flow_t *)realloc(link_flows,
               net->num_links * sizeof(link_data_t))))
      {
        fprintf(stderr, "realloc link_flows in warmstart failed\n");
        return(1);
      }
      while (istart < flows_num_links)
      {
        for (i = istart; i < flows_num_links &&
               link_flows[i].init_node == net->links[i].init_node &&
               link_flows[i].term_node == net->links[i].term_node;
             i++)
          /*nothing*/;
        if (i < flows_num_links)
        {
          /* now the net->links[i] entry is a new one not in link_flows 
             so move everything else "up one" and put a zero flow entry for
             the new link in link_flows */
          memmove(&link_flows[i+1], &link_flows[i], 
                  (flows_num_links - i + inserted_count) * sizeof(link_flow_t));
          link_flows[i].init_node = net->links[i].init_node;
          link_flows[i].term_node = net->links[i].term_node;
          link_flows[i].volume = 0;
          link_flows[i].cost = 0;
          ++inserted_count;
        }
        istart = i + 1;
      }
    }
    else if (flows_num_links > net->num_links)
    {
      /* TODO handle link deletion */
      fprintf(stderr, "error: %ld links in net file but %ld in flows input file\n", net->num_links, flows_num_links);
      return(1);
    }
  }
    


  if (!(Va = (long *)malloc((net->num_nodes+2)*sizeof(long))))
  {
    fprintf(stderr, "malloc Va failed\n");
    return(1);
  }
  if (!(Ea = (long *)malloc(net->num_links*sizeof(long))))
  {
    fprintf(stderr, "malloc Ea failed\n");
    return(1);
  }

  if (!(link_volumes = (double *)malloc(net->num_links * sizeof(double))))
  {
    fprintf(stderr, "malloc link_volumes failed\n");
    return(1);
  }
  if (!(link_volumes2 = (double *)malloc(net->num_links * sizeof(double))))
  {
    fprintf(stderr, "malloc link_volumes2 failed\n");
    return(1);
  }

  if (!(link_costs = (double *)malloc(net->num_links *sizeof(double))))
  {
    fprintf(stderr, "malloc link_costs failed\n");
    return(1);
  }

  if (!(dist = (double *)malloc(num_threads*(net->num_nodes+1) * sizeof(double))))
  {
    fprintf(stderr, "malloc dist failed\n");
    return(1);
  }

  if (!(queue_next = (long *)malloc(num_threads*(net->num_nodes+1)*sizeof(long))))
  {
    fprintf(stderr, "malloc queue failed\n");
    return(1);
  }

  if (!(predlink = (long *)malloc(num_threads*(net->num_nodes+1)*sizeof(long))))
  {
    fprintf(stderr, "malloc pred failed\n");
    return(1);
  }

  /* fix up for 0 capacity otherwise get div. by zero and other problems */
  min_nonzero_capacity = FLOATINF;
  for (i = 0; i < net->num_links; i++)
  {
    if (net->links[i].capacity > 0.0 && net->links[i].capacity < min_nonzero_capacity)
      min_nonzero_capacity = net->links[i].capacity;
  }
  for (i = 0; i < net->num_links; i++)
  {
    if (net->links[i].capacity <= 0.0)
    {
      fprintf(stderr, "warning: zero capacity for link from %ld to %ld, set to min nonzero capacity %f\n", net->links[i].init_node, net->links[i].term_node,
              min_nonzero_capacity);
      net->links[i].capacity = min_nonzero_capacity;
    }
  }

	  
  adjlist_to_packed_arrays(net->links, net->num_nodes+1, net->num_links, Va, Ea, link_costs);

  total_demand = demand_sum(demands, net->num_zones+1);
	    

  if (warm_start_mode)
  {
    /* we sorted the link flow entries so can copy directly to volumes */
    for (i = 0; i < net->num_links; i++) // TODO: parallelize
      link_volumes[i] = link_flows[i].volume;
    /* compute costs based on these volumes, they SHOULD match the ones
       in the file */
    update_link_costs(net, link_volumes, link_costs); // TODO: paralleize
    for (i = 0; i < net->num_links; i++)
    {
      if (fabs(link_costs[i] - link_flows[i].cost) > EPS)
        fprintf(stderr, "warning: mismatched costs (delta = %g) on link %d from %ld to %ld\n", fabs(link_costs[i] - link_flows[i].cost),
                i, link_flows[i].init_node, link_flows[i].term_node);
    }
  }
  else
  {
    for (i = 0; i < net->num_links; i++) // TODO: parallelize
    {
      link_volumes[i] = 0;
      link_volumes2[i] = 0;
    }
  }


  if (!warm_start_mode)
  {
    relative_gap = FLOATINF;
    iter = 0;
    fw_assign(net, demands, Va, Ea, link_costs, link_volumes,
              dist, queue_next, predlink);
    update_link_costs(net, link_volumes, link_costs); // TODO: parallelize
    best_lowerbound = -FLOATINF;
  }

  while ( (target_iterations <= 0 || iter < target_iterations) &&
          (target_relgap <= 0 || relative_gap > target_relgap) &&
          !finished )
  {
    if (time_each_iteration)
      gettimeofday(&start_timeval, NULL);

    fw_assign(net, demands, Va, Ea, link_costs, link_volumes2,
              dist, queue_next, predlink);
    for (i = 0; i < net->num_links; i++) // TODO: parallelize
      link_volumes2[i] = link_volumes2[i] - link_volumes[i]; /* delta volumes */
    lambda = line_search(net, link_volumes, link_volumes2);

    gap = -link_directional_derivative(net, link_volumes, link_volumes2, 0);
    objective_value = links_objective_function(net, link_volumes);
    lowerbound = objective_value - gap;
    if (lowerbound > best_lowerbound)
      best_lowerbound = lowerbound;
    if (best_lowerbound != 0)
      relative_gap = (objective_value - best_lowerbound) / fabs(best_lowerbound);
    // TODO: parallelize total_link_cost() (sum reduction)
    average_excess_cost =  -total_link_cost(net, link_volumes2, link_costs) / 
      total_demand;
    link_cost_total = total_link_cost(net, link_volumes, link_costs);

    for (i = 0; i < net->num_links; i++) // TODO: parallelize
      link_volumes[i] += lambda * link_volumes2[i];
    update_link_costs(net, link_volumes, link_costs); // TODO: paralleize

    if (time_each_iteration)
    {
      gettimeofday(&end_timeval, NULL);
      timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
      etime = 1000*elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
    }

    iter++;

    if (time_each_iteration)
      fprintf(stderr, "iter = %d objective value = %.15f relative gap = %.15f average exccess cost = %.15f total cost = %.15f (%d ms)\n", iter,
              objective_value, relative_gap,
              average_excess_cost, link_cost_total, etime);

    switch (signal_status)
    {
      case SIGNAL_STATUS_NONE:  /* no signal, just continue */
        break;

      case SIGNAL_STATUS_WRITEMSG: /* write statu smsg and continue */ 
        fprintf(stderr, "iter = %d objective value = %.15f relative gap = %.15f average exccess cost = %.15f total cost = %.15f\n", iter,
                objective_value, relative_gap,
                average_excess_cost, link_cost_total);
        break;

      case SIGNAL_STATUS_CHECKPOINT: /* writemsg and flows file and continue */
        fprintf(stderr, "iter = %d objective value = %.15f relative gap = %.15f average exccess cost = %.15f total cost = %.15f\n", iter,
                objective_value, relative_gap,
                average_excess_cost, link_cost_total);
        fprintf(stderr, "checkpointing to flows output file...");
        now = time(NULL);
        rewind(flow_fp);
        fprintf(flow_fp, "~ checkpoint written at %s\n", ctime(&now));
        fprintf(flow_fp,  
                "~ iterations = %d objvalue = %.15f relgap = %.15f\n", 
                iter, objective_value, relative_gap);
        total_cost = total_link_cost(net, link_volumes, link_costs);
        fprintf(flow_fp, "~ totalcost = %.15f\n", total_cost);
        print_flow(flow_fp, net, link_volumes, link_costs);
        fflush(flow_fp);
        fprintf(stderr, "done\n");
        break;

      case SIGNAL_STATUS_TERMINATE:  /* write output and terminate cleanly */
        finished = 1;
        break;

      default:
        fprintf(stderr, "ERROR: unknown signal status %d\n", signal_status);
        break;
    }
    signal_status = SIGNAL_STATUS_NONE;
  }

  total_cost = total_link_cost(net, link_volumes, link_costs);

  rewind(flow_fp);
  fprintf(flow_fp,  
          "~ iterations = %d objvalue = %.15f relgap = %.15f\n", 
          iter, objective_value, relative_gap);
  fprintf(flow_fp, "~ totalcost = %.15f\n", total_cost);
  print_flow(flow_fp, net, link_volumes, link_costs);
  fflush(flow_fp);

  /* cleanup and return */
  free(link_costs);
  free(link_volumes);
  free(link_volumes2);
  free(Ea);
  free(Va);
  for (i = 0; i < net->num_zones+1; i++)
    free(demands[i]);
  free(demands);
  free(queue_next);
  free(predlink);
  free(dist);
  
  *total_cost_at_end = total_cost;
  *objvalue_at_end = objective_value;
  *iterations_at_end = iter;
  *relgap_at_end = relative_gap;

  return 0;

}
