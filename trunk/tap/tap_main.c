/*****************************************************************************
 * 
 * File:    tap_main.c
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_main.c 787 2011-10-05 23:41:39Z astivala $
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
 *                                 [-t tripmodfilename] [-v]
 *                                 netfilename  demandfilename
 *
 *     netfilename    is name of the net file defining node and links
 *     demandfilename is name of Origin-Destination demand file
 *     -n threads       : run with that many worker threads
 *     -w flowsfilename : warmstart from flows stored in flowsfilename
 *                        which is the output link flows file format
 *     -i iterations    : terminate after this many iterations
 *     -r relgap        : terminate after relative gap is <= relgap.
 *    -t tripmodfilename: parse modifications to the O-D demand data form
 *                        tripmodfilename and run with the modified trip data
 *     -v               : verbose: output objective value, relgap, time, etc. each iteration
 *                        to stderr
 *                        
 *     Output is link flows on stdout.
 *     
 *   Example usage:
 *   
 *   tap_frankwolfe_pthread SiouxFalls_net.txt SiouxFalls_trips.txt > SioxFalls_flow.txt
 *
 * Three signals are handled:
 *
 *     SIGUSR1 (signal 10) - write status (current values, iteration etc.)
 *                           to stderr and continue.
 *
 *     SIGUSR2 (signal 12) - write status (current values, iteration etc.)
 *                           to stderr and also write flow output file.
 *                           The -w option can be used to continue later
 *                           as a kind of checkpoint/restart.
 *
 *     SIGHUP  (signal 1)  - terminate cleanly at end of current iteration,
 *                           writing flow output file etc. just as if
 *                           the iteration or relative gap limit was reached.
 *                           The -w option can be used to continue later.
 *   
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#include <pthread.h>
#include <unistd.h>
#include "utils.h"
#include "parsetapfiles.h"
#include "tap_frankwolfe_pthread.h"
#include "tap_functions.h"
   
/*****************************************************************************
 *
 * constants
 *
 ****************************************************************************/

const int DEFAULT_ITERATION_LIMIT = 10000;


/*****************************************************************************
 *
 * Main
 *
 ****************************************************************************/

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s  [-n numthreads] [-w warmstart_flows_filename] [-t tripmods_filename] netfilename demandfilename\n"
          "  -n numthreads : number of threads to use (default %d)\n"
          "  -w warmsetart_flows_filename : warmstart from previous solution\n"
          "  -i iterations : terminate after this many iterations (default %d)\n"
          "  -r relgap     : terminate when relative gap <= relgap\n"
          "  -t tripmodsfilename : modifiy O-D demand data\n"
          "  -v : verbose - write status each iteration to stderr\n"
          , progname, num_cores, DEFAULT_ITERATION_LIMIT );
  exit(1);
}

int main(int argc, char *argv[])
{
  char *net_filename, *demand_filename;
  FILE *net_fp, *demand_fp;
  int c;
  FILE *flows_input_fp = NULL;
  char *flows_input_filename = NULL;
  int target_iterations = 0;
  double target_relgap = 0;
  int warm_start_mode = 0;
  double total_cost;
  net_data_t net;
  int iterations;
  double relgap,objvalue;
  char *trip_mods_filename = NULL;
  FILE *trip_mods_fp = NULL;
  int trip_mods_mode = 0;
  int num_trip_mods = 0;
  trip_mod_t *trip_mods = NULL;
  demand_data_t **demands = NULL;
  long num_zones = 0;
  int k;
  int verbose = 0;

  num_cores = get_num_cores();
  num_threads = num_cores;

  while ((c = getopt(argc, argv, "n:w:i:r:t:v")) != -1)
  {
    switch (c)
    {
      case 'n':  /* number of threads */
        if (atoi(optarg) < 1)
        {
          fprintf(stderr, "number of threads must be >= 1\n");
          usage(argv[0]);
        }
        else if (atoi(optarg ) > MAX_NUM_THREADS)
        {
          fprintf(stderr, "maximum number of threads is %d\n", MAX_NUM_THREADS);
          usage(argv[0]);
        }
        else
          num_threads = atoi(optarg);
        break;

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

      case 't': /* modifi origin-dest demand data not network */
        trip_mods_filename = optarg;
        trip_mods_mode = 1;
        break;

      case 'v': /* verbose: write status to stderr each iteration */
        verbose = 1;
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

  if (trip_mods_mode)
  {
    if (!(trip_mods_fp = fopen(trip_mods_filename, "r")))
    {
      fprintf(stderr, "error opening trips modifications file %s: %s\n",
              trip_mods_filename, strerror(errno));
      exit(1);
    }
    if (parse_trip_mod_file(trip_mods_fp, &trip_mods, &num_trip_mods) != 0)
    {
      fprintf(stderr, "parse trip mods file %s failed\n", trip_mods_filename);
      exit(1);
    }
    fclose(trip_mods_fp);
  }

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
  }

  if (parse_net_file(net_fp, &net) != 0)
  {
    fprintf(stderr, "error parsing net file\n");
    exit(1);
  }
  fclose(net_fp);

  /* sort link data so matches saved flow data if warmstart used */
  qsort(net.links, net.num_links, sizeof(link_data_t), link_data_compar);

  if (parse_trips_file(demand_fp, &demands, &num_zones) != 0)
  {
    fprintf(stderr, "error parsing trips file\n");
    return(1);
  }
  fclose(demand_fp);
  if (num_zones != net.num_zones)
  {
    fprintf(stderr, "warning: %ld zones in net data but %ld in trips file\n",
            net.num_zones, num_zones);
  }

  if (trip_mods_mode)
  {
    for (k = 0; k < num_trip_mods; k++)
      apply_trip_modification(demands, &trip_mods[k], num_zones);
  }

  if (tap_frankwolfe(&net, demands, warm_start_mode, 
                     flows_input_fp,
                     target_iterations, target_relgap, stdout,
                     &total_cost, &objvalue, &iterations, &relgap, verbose) != 0)
  {
    fprintf(stderr, "failed\n");
    exit(1);
  }
  fprintf(stderr, "iterations = %d objective value = %.15f relative gap = %.15f total cost = %.15f\n", iterations, objvalue, relgap, total_cost);

  exit(0);
}

