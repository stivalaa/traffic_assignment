/*****************************************************************************
 * 
 * File:    tap_cuda_main.c
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_cuda_main.cu 702 2011-09-14 23:26:53Z astivala $
 *
 * Traffic assignment by Frank-Wolfe algorithm.
 * This version runs the shortest path computations (d'Esopo-Pape 
 * algorithm) in parallel ona GPU using CUDA
 *
 * There are many descriptions of this algorithm, a simple one for
 * example is in Sheffi, Y. 1985 "Urban Transportation Networks:
 * Equilibribum Analysis with Mathematical Programming Methods"
 * Prentice-Hall (out of print, available online from
 * http://web.mit.edu/sheffi/www/urbanTransportation.html
 *
 *   Usage: tap_frankwolfe_cuda  [-w flowsfilename] 
 *                                 [-i iter] [-r relgap] 
 *                                 [-t tripmodfilename]
 *                                 netfilename  demandfilename
 *
 *     netfilename    is name of the net file defining node and links
 *     demandfilename is name of Origin-Destination demand file
 *     -w flowsfilename : warmstart from flows stored in flowsfilename
 *                        which is the output link flows file format
 *     -i iterations    : terminate after this many iterations
 *     -r relgap        : terminate after relative gap is <= relgap.
 *    -t tripmodfilename: parse modifications to the O-D demand data form
 *                        tripmodfilename and run with the modified trip data
 *
 *
 *     Output is link flows on stdout.
 *     
 *   Example usage:
 *   
 *   tap_frankwolfe_cuda SiouxFalls_net.txt SiouxFalls_trips.txt > SioxFalls_flow.txt
 *
 *
 ****************************************************************************/
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include "utils.h"
#include "parsetapfiles.h"
#include "tap_frankwolfe_cuda.h"
#include "tap_functions.h"
#include "volume_update_host.h"
 
#define VERBOSE
  
/* FIXME -  should either fix the start_nodes vector code to actually allow arbitrary start nodes or remove it and hardcode to have start nodes 1..num_zones  which it currently assumes (with hack to add 1 s for not 0.. num_zones-1) */

/*****************************************************************************
 *
 * constants
 *
 ****************************************************************************/

const int DEFAULT_ITERATION_LIMIT = 10000;


/*****************************************************************************
 *
 * globals
 *
 ****************************************************************************/
int using_fermi_architecture = 0; // using a Fermi architecture GPU


/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/


/*
 * choose_gpu() - choose the GPU to use
 *
 * Parameters: 
 *      major_req - required (minimum) major compute capability
 *      minor_req - required (minimum) minor compute capabililty
 * Return value: None
 *
 * Sets the global using_fermi_architecture and calls cudaSetDevice()
 * to set the GPU to use.
 */
static void choose_gpu(int major_req, int minor_req)
{
/*
    long devnum = cutGetMaxGflopsDeviceId();
    fprintf(stderr, "using max gflops device %d: ", devnum);
*/
    /* If there is a compute capability 2 device ("Fermi"
       architecture) (or higher) then use that
    */

#if defined(__DEVICE_EMULATION__)
    fprintf(stderr, "device emulation on\n");
    cudaSetDevice( 0 );
    return;
#endif

    int devnum, deviceCount, gflops,max_gflops=0, sel_devnum=-1;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
    {
      fprintf(stderr, "There is no device supporting CUDA.\n");
      exit(1);
    }
    for (devnum = 0; devnum < deviceCount; devnum++)
    {  
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, devnum);
      if (devnum == 0 && 
          deviceProp.major == 9999 && deviceProp.minor == 9999)
      {
        fprintf(stderr, "There is no device supporting CUDA.\n");
        exit(1);
      }

      fprintf(stdout, "found %d CUDA devices\n", deviceCount);

      if (deviceProp.major < major_req ||
          deviceProp.major == major_req && deviceProp.minor < minor_req)
      {
        fprintf(stderr,
             "cannot use device %d compute capability %d.%d (require %d.%d)\n",
             devnum, deviceProp.major, deviceProp.minor, major_req, minor_req);
        
        continue;
      }

      if (deviceProp.major >= 2)
      {
        fprintf(stdout,
          "found Fermi architecture (compute capability %d.%d) device %d: %s\n",
                deviceProp.major, deviceProp.minor, devnum, deviceProp.name);
        sel_devnum = devnum;
        using_fermi_architecture = 1;
        break;
      }
      else
      {
        gflops = deviceProp.multiProcessorCount * deviceProp.clockRate;
        fprintf(stdout, "device %d: %s\n", devnum,
                deviceProp.name);
        if (gflops > max_gflops)
        {
          max_gflops = gflops;
          sel_devnum = devnum;
          using_fermi_architecture = 0;
        }
      }
    }
    if (sel_devnum < 0)
    {
      fprintf(stderr,
              "there are no CUDA devices of required compute capability\n");
      exit(1);
    }

    fprintf(stdout, "using device %d: ", sel_devnum);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, sel_devnum);
    fprintf(stdout, "%s\n", deviceProp.name);
    cudaSetDevice( sel_devnum );
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
  long *predlink;
  int target_iterations = 0;
  double target_relgap = 0;

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

  /* Initialize GPU */
  // need at least compute capability 1.3 for double precision f.p.
  choose_gpu(1,3);  /* calls cudaSetDevice() */

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


  if (!(dist = (double *)malloc((net.num_nodes+1) * (1+net.num_zones )* sizeof(double))))
  {
    fprintf(stderr, "malloc dist failed\n");
    exit(1);
  }

  if (!(predlink = (long *)malloc((net.num_nodes+1)*(1+net.num_zones)*sizeof(long))))
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

#ifdef USE_CUDA_VOLUME_UPDATE
  link_volume_data_setup(net.links, net.num_links, net.num_zones, demands);
#endif

  if (!warm_start_mode)
  {
    relative_gap = FLOATINF;
    iter = 0;
    fw_assign(&net, demands, Va, Ea, link_costs, link_volumes, dist,  predlink);

    update_link_costs(&net, link_volumes, link_costs);
    best_lowerbound = -FLOATINF;
  }

  while ( (target_iterations <= 0 || iter < target_iterations) &&
          (target_relgap <= 0 || relative_gap > target_relgap) )
  {
    gettimeofday(&start_timeval, NULL);
    fw_assign(&net, demands, Va, Ea, link_costs, link_volumes2,dist,
              predlink);
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
  free(predlink);
#ifdef USE_CUDA_VOLUME_UPDATE
  link_volume_data_cleanup();
#endif
  cudaThreadExit();
  exit(0);
}

