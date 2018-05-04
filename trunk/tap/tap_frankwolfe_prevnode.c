/*****************************************************************************
 * 
 * File:    tap_frankwolfe.c
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_frankwolfe_prevnode.c 479 2011-07-13 07:14:16Z astivala $
 *
 * Traffic assignment by Frank-Wolfe algorithm.
 *
 * There are many descriptions of this algorithm, a simple one for
 * example is in Sheffi, Y. 1985 "Urban Transportation Networks:
 * Equilibribum Analysis with Mathematical Programming Methods"
 * Prentice-Hall (out of print, available online from
 * http://web.mit.edu/sheffi/www/urbanTransportation.html
 *
 *   Usage: tap_frankwolfe [-w flowsfilename] [-i iter] [-r relgap]
 *                          netfilename  demandfilename
 *
 *     netfilename    is name of the net file defining node and links
 *     demandfilename is name of Origin-Destination demand file
 *     -w flowsfilename : warmstart from flows stored in flowsfilename
 *                        which is the output link flows file format
 *     -i iterations    : terminate after this many iterations
 *     -r relgap        : terminate after relative gap is <= relgap.
 *
 *     Output is link flows on stdout.
 *     
 *   Example usage:
 *   
 *   tap_frankwolfe SiouxFalls_net.txt SiouxFalls_trips.txt > SioxFalls_flow.txt
 *   
 * TODO merge in [warm start] added link capability frmo tap_frankwolfe_pthread.c
 *      and in general work out what to do with this version- mostly
 *      duplicated code, but might want to retain it since there is
 *      (sometimes significant) overhead using pthreads at all for even 1
 *      thread versus this unthreaded version.
 *
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include "parsetapfiles.h"
#include "sssp_pape.h"
#include "tap_functions.h"
#include "utils.h"


#define VERBOSE

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
   are not numbered externally.  Hence we are always adding 1 to
   num_nodes in allocating vectors.
*/
   
   
/*****************************************************************************
 *
 * constants
 *
 ****************************************************************************/

const int DEFAULT_ITERATION_LIMIT = 10000;



/*****************************************************************************
 *
 * Local functions
 *
 ****************************************************************************/

/*
 *  fw_assign() - flow assignment in Frank-Wolfe
 *
 *   Assign flow from O-D demands to links according to current shortest
 *   paths (lowest costs) on netgraph
 *
 * Parameters:
 *    net - net structure parsed from net file
 *    demands  - demands array parsed from trips file
 *                for each origin demands[origin] is an array of demands structs
 *                terminated by one with 0  dest (origins and dests start at 1)
 *    Va, Ea, Wa - network in packed adjancey list format (Wa is current costs)
 *    link_volumes - (IN/OUT) - volume on each link
 *   dist (WORK) - vector of num_nodes doubles for distance from origin of each
 *                 for sssp_pape()
 *   queue_next (WORK) - vector of num_nodes longs for queue of nodes
 *                 for sssp_pape()
 *   predlink (WORK) -vector of num_nodes longs for predecessor link each node
 *   prednode - (WORK) vector of num_nodes longs for predecessor of each node
 *   old_prednodes - (WORK) vectors of num_nodes longs for
 *                    predecessor of each node, one for each origin
 *
 *
 *  Return value;
 *    None.
 */
static void fw_assign(net_data_t *net, demand_data_t *demands[],
                      long Va[], long Ea[], double Wa[], 
                      double link_volumes[],
                      double dist[], long queue_next[],
                      long predlink[],
                      long prednode[], long old_prednodes[])
{
  long orig, dest;
  int i;
  double routeflow;
  long v,k;
  int fail;

#ifdef VERIFY_SSSP
  long predlink2[100000];
  double dist2[100000]; 
  double pathcost;
  int j;
#endif

  for (k = 0; k < net->num_links; k++)
    link_volumes[k] = 0;

  for (orig = 1; orig < net->num_zones+1; orig++)
  {
    /* TODO the prevnodefirst heuristic was not much better than simpler
       pape_lll version, should just revert to that (see notes.txt) */
    sssp_prevnodefirst_lll(Va, Ea, Wa, net->num_nodes+1, net->num_links, 
                       orig, net->first_thru_node, 
                       &old_prednodes[orig * (net->num_nodes+1)],
                       predlink, prednode,
                       dist, queue_next);
    for (i = 0; i < net->num_nodes+1; i++) /*copy prednode to old for next iter*/
      old_prednodes[orig * (net->num_nodes+1) + i] = prednode[i];

#ifdef VERIFY_SSSP
    sssp_pape(Va, Ea, Wa, net->num_nodes+1, net->num_links, 
              orig, net->first_thru_node, predlink2, dist2, queue_next);
    for (j = 0 ; j <net->num_nodes+1; j++)
    {
      if (fabs(dist[j] - dist2[j]) > 1e-04)
        fprintf(stderr, "ERROR at j = %d, dist[j] = %f but dist2[j] = %f\n",j,dist[j],dist2[j]);
    }
    for (j = 0; (dest = demands[orig][j].dest) != 0; j++)
    {
      if (dest == orig) continue;
      pathcost = 0;
      v = dest;
      while (v != orig)
      {
        k = predlink[v];
        assert(k >= 0);
        assert(net->links[k].term_node == v);
        if (predlink2[v] != k)
           fprintf(stderr, "WARN path differs for orig = %d dest = %d at v = %d\n", j, dest, v);
        pathcost += Wa[k];
        v = net->links[k].init_node;
      }
      if (fabs(pathcost - dist[dest] > 1e-04))
        fprintf(stderr, "ERROR at orig = %d dest = %d, dist[dest] = %f but pathcost = %f\n", j, dest, dist[dest], pathcost);
    }
#endif

    for (i = 0; (dest = demands[orig][i].dest) != 0; i++)
    {
      routeflow = demands[orig][i].demand;
      if (orig == dest || routeflow <= 0)/*should never be -ve; stops warning*/
        continue;
      v = dest;
      fail = 0;
      while (v != orig)
      {
        k = predlink[v];
        if (k < 0)
        {
          /* this shouldn't happen, it means there is a nonzero demand
             between two nodes but no path between them  - 
             probably due to bad input data */
          fprintf(stderr, "ERROR: no path from %ld to %ld\n", orig, dest);
          exit(1);
          break;
        }
        assert(k >= 0);
        assert(net->links[k].term_node == v);
        assert(prednode[v] == net->links[k].init_node);
        link_volumes[k] += routeflow;
        v = net->links[k].init_node;
      }
    }
  }
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
 * Main
 *
 ****************************************************************************/

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s  [-w warmstart_flows_filename] netfilename demandfilename\n"
          "  -w warmsetart_flows_filename : warmstart from previous solution\n"
          "  -i iterations : terminate after this many iterations (default %d)\n"
          "  -r relgap     : terminate when relative gap <= relgap\n"
          , progname,  DEFAULT_ITERATION_LIMIT );
  exit(1);
}

int main(int argc, char *argv[])
{
  char *net_filename, *demand_filename;
  FILE *net_fp, *demand_fp;
  net_data_t net;
  demand_data_t **demands = NULL;
  long num_zones = 0;
  long *Va,*Ea;  /* vertex and edge arrays for packed adjlist format of graph */
  int i,j;
  double *link_volumes, *link_costs, *link_volumes2;
  int iter;
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int etime;
  double min_nonzero_capacity;
  double lambda;
  double total_demand = 0;
  double gap, relative_gap, lowerbound, best_lowerbound, objective_value;
  double link_cost_total, average_excess_cost;
  int c;
  FILE *flows_input_fp = NULL;
  char *flows_input_filename = NULL;
  long flows_num_links = 0;
  link_flow_t *link_flows = NULL;
  int warm_start_mode = 0;
  double *dist;
  long *queue_next, *predlink;
  int target_iterations = 0;
  double target_relgap = 0;
  long *prednode, *old_prednodes;

  while ((c = getopt(argc, argv, "w:i:r:")) != -1)
  {
    switch (c)
    {
      case 'w':  /* warmstart file */
        flows_input_filename = optarg;
        warm_start_mode = 1;
        break;
        
      case 'i': /* iteration limit */
        if (atoi(optarg) < 1)
        {
          fprintf(stderr, "iteration target must be >= 1\n");
          usage(argv[0]);
        }
        target_iterations = atoi(optarg);
        break;

      case 'r': /* relative gap target */
        if (atof(optarg) <= 0)
        {
          fprintf(stderr, "relative gap must be positive\n");
          usage(argv[0]);
        }
        target_relgap = atof(optarg);
        break;

      default:
        usage(argv[0]);
        break;
    }
  }
  
  if (argc - optind != 2)
    usage(argv[0]);

  net_filename = argv[optind];
  demand_filename = argv[optind+1];

  if (target_relgap <= 0 && target_iterations <= 0)
    target_iterations = DEFAULT_ITERATION_LIMIT;

  if (!(net_fp = fopen(net_filename, "r")))
  {
    fprintf(stderr, "error opening net file %s: %s\n", 
            net_filename, strerror(errno));
    exit(1);
  }
  
  if (!(demand_fp = fopen(demand_filename, "r")))
  {
    fprintf(stderr, "error opening trips file %s: %s\n",
            demand_filename, strerror(errno));
    exit(1);
  }

  if (parse_net_file(net_fp, &net) != 0)
  {
    fprintf(stderr, "error parsing net file %s\n", net_filename);
    exit(1);
  }
  fclose(net_fp);

  if (parse_trips_file(demand_fp, &demands, &num_zones) != 0)
  {
    fprintf(stderr, "error parsing trips file %s\n", demand_filename);
    exit(1);
  }
  fclose(demand_fp);

  if (warm_start_mode)
  {
    /* we ahve the 'warm start' option to parse link volumes (and costs)
       from a flows file (output of this program also) */
    if (!(flows_input_fp = fopen(flows_input_filename, "r")))
    {
      fprintf(stderr, "error opening flows input file %s: %s\n",
              flows_input_filename, strerror(errno));
      exit(1);
    }
    if (parse_flows_file(flows_input_fp, &link_flows, &flows_num_links) != 0)
    {
      fprintf(stderr, "error parsing flows input file %s\n", flows_input_filename);
      exit(1);
    }
    fclose(flows_input_fp);
    if (flows_num_links != net.num_links)
    {
      fprintf(stderr, "error: %ld links in net file but %ld in flows input file\n", net.num_links, flows_num_links);
      exit(1);
    }
    /* sort by 'from' node ascending and within that by 'to' node ascending */
    qsort(link_flows,flows_num_links,sizeof(link_flow_t),linkflow_entry_compar);
  }
    

  if (num_zones != net.num_zones)
  {
    fprintf(stderr, "warning: %ld zones in net data but %ld in trips file\n",
            net.num_zones, num_zones);
  }

  if (!(Va = (long *)malloc((net.num_nodes+2)*sizeof(long))))
  {
    fprintf(stderr, "malloc Va failed\n");
    exit(1);
  }
  if (!(Ea = (long *)malloc(net.num_links*sizeof(long))))
  {
    fprintf(stderr, "malloc Ea failed\n");
    exit(1);
  }

  if (!(link_volumes = (double *)malloc(net.num_links * sizeof(double))))
  {
    fprintf(stderr, "malloc link_volumes failed\n");
    exit(1);
  }
  if (!(link_volumes2 = (double *)malloc(net.num_links * sizeof(double))))
  {
    fprintf(stderr, "malloc link_volumes2 failed\n");
    exit(1);
  }

  if (!(link_costs = (double *)malloc(net.num_links *sizeof(double))))
  {
    fprintf(stderr, "malloc link_costs failed\n");
    exit(1);
  }


  if (!(dist = (double *)malloc((net.num_nodes+1) * sizeof(double))))
  {
    fprintf(stderr, "malloc dist failed\n");
    exit(1);
  }

  if (!(queue_next = (long *)malloc((net.num_nodes+1) *sizeof(long))))
  {
    fprintf(stderr, "malloc queue failed\n");
    exit(1);
  }

  if (!(predlink = (long *)malloc((net.num_nodes+1)*sizeof(long))))
  {
    fprintf(stderr, "malloc pred failed\n");
    exit(1);
  }

  if (!(prednode = (long *)malloc((net.num_nodes+1)*sizeof(long))))
  {
    fprintf(stderr, "malloc prednode failed\n");
    exit(1);
  }
  if (!(old_prednodes = (long *)malloc(
          (net.num_zones+1)*(net.num_nodes+1)*sizeof(long))))
  {
    fprintf(stderr, "malloc old_prednodes failed\n");
    exit(1);
  }

  /* fix up for 0 capacity otherwise get div. by zero and other problems */
  min_nonzero_capacity = FLOATINF;
  for (i = 0; i < net.num_links; i++)
  {
    if (net.links[i].capacity > 0.0 && net.links[i].capacity < min_nonzero_capacity)
      min_nonzero_capacity = net.links[i].capacity;
  }
  for (i = 0; i < net.num_links; i++)
  {
    if (net.links[i].capacity <= 0.0)
    {
      fprintf(stderr, "warning: zero capacity for link from %ld to %ld, set to min nonzero capacity %f\n", net.links[i].init_node, net.links[i].term_node,
        min_nonzero_capacity);
      net.links[i].capacity = min_nonzero_capacity;
    }
  }

  
  adjlist_to_packed_arrays(net.links, net.num_nodes+1, net.num_links, Va, Ea, link_costs);

  total_demand = demand_sum(demands, net.num_zones+1);
    

  if (warm_start_mode)
  {
    /* we sorted the link flow entries so can copy directly to volumes */
    for (i = 0; i < net.num_links; i++)
      link_volumes[i] = link_flows[i].volume;
    /* compute costs based on these volumes, they SHOULD match the ones
       in the file */
    update_link_costs(&net, link_volumes, link_costs);
    for (i = 0; i < net.num_links; i++)
    {
      if (fabs(link_costs[i] - link_flows[i].cost) > EPS)
        fprintf(stderr, "warning: mismatched costs (delta = %g) on link %d from %ld to %ld\n", fabs(link_costs[i] - link_flows[i].cost),
                i, link_flows[i].init_node, link_flows[i].term_node);
    }
  }
  else
  {
    for (i = 0; i < net.num_links; i++)
    {
      link_volumes[i] = 0;
      link_volumes2[i] = 0;
    }
  }


  for (i = 0; i < net.num_nodes + 1; i++)
  {
    prednode[i] = -1;    /* INVALID */
  }
  for (i = 0; i < net.num_zones + 1; i++)
  {
    for (j = 0; j < net.num_nodes + 1; j++)
    {
      old_prednodes[i*(net.num_zones+1)+j] = -1;   /* INVALID */
    }
  }

    
  if (!warm_start_mode)
  {
    relative_gap = FLOATINF;
    iter = 0;
    fw_assign(&net, demands, Va, Ea, link_costs, link_volumes, dist, queue_next,
              predlink, prednode, old_prednodes);
    update_link_costs(&net, link_volumes, link_costs);
    best_lowerbound = -FLOATINF;
  }

  while ( (target_iterations <= 0 || iter < target_iterations) &&
          (target_relgap <= 0 || relative_gap > target_relgap) )
  {
    gettimeofday(&start_timeval, NULL);
    fw_assign(&net, demands, Va, Ea, link_costs, link_volumes2,dist, queue_next,
              predlink, prednode, old_prednodes);
    for (i = 0; i < net.num_links; i++)
      link_volumes2[i] = link_volumes2[i] - link_volumes[i]; /* delta volumes */
    lambda = line_search(&net, link_volumes, link_volumes2);

    gap = -link_directional_derivative(&net, link_volumes, link_volumes2, 0);
    objective_value = links_objective_function(&net, link_volumes);
    lowerbound = objective_value - gap;
    if (lowerbound > best_lowerbound)
      best_lowerbound = lowerbound;
    if (best_lowerbound != 0)
      relative_gap = (objective_value - best_lowerbound) / fabs(best_lowerbound);
    average_excess_cost =  -total_link_cost(&net, link_volumes2, link_costs) / 
      total_demand;
    link_cost_total = total_link_cost(&net, link_volumes, link_costs);

    for (i = 0; i < net.num_links; i++)
      link_volumes[i] += lambda * link_volumes2[i];
    update_link_costs(&net, link_volumes, link_costs);

    gettimeofday(&end_timeval, NULL);
    timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
    etime = 1000*elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
    iter++;
#ifdef VERBOSE
    fprintf(stderr, "iter = %d objective value = %.15f relative gap = %.15f average exccess cost = %.15f total cost = %.15f (%d ms)\n", iter,
            objective_value, relative_gap,
            average_excess_cost, link_cost_total, etime);
#endif
  }

  fprintf(stderr, "total cost = %.15f\n", 
          total_link_cost(&net, link_volumes, link_costs));
  print_flow(stdout, &net, link_volumes, link_costs);

  /* cleanup and exit */
  free(link_costs);
  free(link_volumes);
  free(link_volumes2);
  free(Ea);
  free(Va);
  free(net.links);
  for (i = 0; i < num_zones+1; i++)
    free(demands[i]);
  free(demands);
  free(dist);
  free(queue_next);
  free(predlink);
  free(prednode);
  free(old_prednodes);
  exit(0);
}

