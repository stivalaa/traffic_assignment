/*****************************************************************************
 * 
 * File:    tap_greedy.c
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_greedy.c 245 2011-04-21 00:17:45Z astivala $
 *
 * Traffic assignment by greedy algorithm.
 *
 *   Usage: tap_greedy netfilename  demandfilename
 *
 *     netfilename    is name of the net file defining node and links
 *     demandfilename is name of Origin-Destination demand file
 *
 *     Output is link flows on stdout.
 *     
 *   Example usage:
 *   
 *   tap_greedy SiouxFalls_net.txt SiouxFalls_trips.txt > SioxFalls_flow.txt
 *   
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "parsetapfiles.h"
#include "sssp_pape.h"
#include "tap_functions.h"

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
 * Local functions
 *
 ****************************************************************************/

/*
 *  greedy_assign() - greedy flow assignment
 *
 *   Assign flow from O-D demands to links according to current shortest
 *   paths (lowest costs) on netgraph
 *   Assigns flow to path only up to the amount of flow that
 *   takes the first link to do so up to the next step in cost function,
 *   subtracting that assigned flow from the O-D demand using that path
 *
 * Parameters:
 *    net - net structure parsed from net file
 *    demands (IN/OUT) - demands array parsed from trips file
 *                for each origin demands[origin] is an array of demands structs
 *                terminated by one with 0  dest (origins and dests start at 1)
 *    Va, Ea, Wa - network in packed adjancey list format (Wa is current costs)
 *    link_volumes - (IN/OUT) - volume on each link
 *   dist (WORK) - vector of num_nodes doubles for distance from origin of each
 *                 for sssp_pape()
 *   queue_next (WORK) - vector of num_nodes longs for queue of nodes
 *                 for sssp_pape()
 *   predlink (WORK) - vector of num_nodes longs for predecessor of each node
 *
 *  Return value;
 *    None.
 */
static void greedy_assign(net_data_t *net, demand_data_t *demands[],
                          long Va[], long Ea[], double Wa[], 
                          double link_volumes[],
                          double dist[], long queue_next[], long predlink[])
{
  long orig, dest;
  int i;
  double routeflow;
  long v,k;
  double addflow;
  int fail;


  for (orig = 1; orig < net->num_zones+1; orig++)
  {
    sssp_pape(Va, Ea, Wa, net->num_nodes+1, net->num_links, 
              orig, net->first_thru_node, predlink, dist, queue_next);

    for (i = 0; (dest = demands[orig][i].dest) != 0; i++)
    {
      routeflow = demands[orig][i].demand;
      if (orig == dest || routeflow <= 0)/*should never be -ve; stops warning*/
        continue;

      /* assign flow to path only up to the amount of flow that takes the
         first link to do so up to the next step in cost function */
      v = dest;
      addflow = routeflow;
      fail = 0;
      while (v != orig) /* this loop finds the flow we can add */
      {
        double flow_to_next_step;
        k = predlink[v];
        if (k < 0)
        {
          int x;/*XXX*/
          /* this shouldn't happen, it means there is a nonzero demand
             between two nodes but no path between them  - 
             probably due to bad input data  -mabe this should actually
             be a fatal error FIXME */
          fprintf(stderr, "warning: no path from %ld to %ld\n", orig, dest);
          /* /\*XXX*\/          fprintf(stderr, "k = %ld v = %ld\n",k,v); */
          /* /\*XXX*\/ for(x = 0; x < net->num_nodes; x++) fprintf(stderr,"%ld ",predlink[x]); printf("\n"); */
          fail = 1;
          break;
        }
        assert(net->links[k].term_node == v);
        flow_to_next_step = distance_from_next_step(net->links, link_volumes, k);
        if (flow_to_next_step < addflow)
          addflow = flow_to_next_step;
        v = net->links[k].init_node;
      }
      if (fail)
        continue; /* next dest */
      v = dest;
      while (v != orig) /* then this loop adds that flow to link volume */
      {
        k = predlink[v];
        assert(k >= 0);
        assert(net->links[k].term_node == v);
        link_volumes[k] += addflow;
        v = net->links[k].init_node;
      }
      demands[orig][i].demand -= addflow; /* and subtracts from demand */
    }
  }
}


/*****************************************************************************
 *
 * Main
 *
 ****************************************************************************/

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s netfilename demandfilename\n", progname);
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
  int i;
  double *link_volumes, *link_costs;
  double unsatisfied_demand;
  int iter;
  long orig,dest;
  struct rusage starttime,endtime;
  int otime;
  double min_nonzero_capacity;
  double prev_unsatisfied_demand;
  double *dist;
  long *queue_next, *predlink;
  
  if (argc != 3)
    usage(argv[0]);

  net_filename = argv[1];
  demand_filename = argv[2];

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

  if (parse_trips_file(demand_fp, &demands, &num_zones) != 0)
  {
    fprintf(stderr, "error parsing trips file %s\n", demand_filename);
    exit(1);
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
  if (!(link_costs = (double *)malloc(net.num_links *sizeof(double))))
  {
    fprintf(stderr, "malloc link_costs failed\n");
    exit(1);
  }
  for (i = 0; i < net.num_links; i++)
    link_volumes[i] = 0;

  adjlist_to_packed_arrays(net.links, net.num_nodes+1, net.num_links, Va, Ea, link_costs);

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

  /* set demand for origin=dest  to 0 otherwise greedy alg. won't terminate */
  for (orig = 1; orig < num_zones+1; orig++)
  {
/*    fprintf(stderr, "orig %ld\n", orig); // XXX */
    for (i = 0; (dest = demands[orig][i].dest) != 0; i++)
    {
/*      fprintf(stderr, "  %ld : %f\n", dest,demands[orig][i].demand);//XXX*/
      if (demands[orig][i].dest == orig)
      {
        fprintf(stderr, 
            "warning: nonzero demand for orig = dest = %ld: set to 0\n", orig);
        demands[orig][i].demand = 0;
      }
    }
  }

  iter = 0;
  unsatisfied_demand = demand_sum(demands, net.num_zones+1);
  fprintf(stderr, "initial unsatifised demand = %f\n", unsatisfied_demand);
  while (unsatisfied_demand > 0)
  {
    getrusage(RUSAGE_SELF, &starttime);
    greedy_assign(&net, demands, Va, Ea, link_costs, link_volumes,
                  dist, queue_next, predlink);
    update_link_costs(&net, link_volumes, link_costs);
    prev_unsatisfied_demand = unsatisfied_demand;
    unsatisfied_demand = demand_sum(demands, net.num_zones+1);
    getrusage(RUSAGE_SELF, &endtime);
    otime = 1000 * endtime.ru_utime.tv_sec + endtime.ru_utime.tv_usec/1000
      + 1000 * endtime.ru_stime.tv_sec + endtime.ru_stime.tv_usec/1000
      - 1000 * starttime.ru_utime.tv_sec - starttime.ru_utime.tv_usec/1000
      - 1000 * starttime.ru_stime.tv_sec - starttime.ru_stime.tv_usec/1000;
    iter++;
    fprintf(stderr, "iter = %d unsatisfied demand = %f (time %d ms)\n", iter,
            unsatisfied_demand, otime);
    if (fabs(unsatisfied_demand - prev_unsatisfied_demand) < EPS)
    {
      /* if we have dodgy input data such as demands between origin 
         and destination with no path, it may be impossible to reduce
         unsatsified demand to 0 but we can detect that it has not
         reduced from the last iteration */
      fprintf(stderr, "WARNING: unsatisfied demand unchanged, terminating with unsatisfied demand = %f\n", unsatisfied_demand);
      break;
    }
  }

  fprintf(stderr, "total cost = %f\n", 
          total_link_cost(&net, link_volumes, link_costs));
  print_flow(stdout, &net, link_volumes, link_costs);

  /* cleanup and exit */
  free(link_costs);
  free(link_volumes);
  free(Ea);
  free(Va);
  free(net.links);
  for (i = 0; i < num_zones+1; i++)
    free(demands[i]);
  free(demands);
  free(dist);
  free(queue_next);
  free(predlink);
  exit(0);
}

