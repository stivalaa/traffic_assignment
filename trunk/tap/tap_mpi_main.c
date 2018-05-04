/*****************************************************************************
 * 
 * File:    tap_mpi_main.c
 * Author:  Alex Stivala
 * Created: May 2011
 *
 * $Id: tap_mpi_main.c 787 2011-10-05 23:41:39Z astivala $
 *
 * Traffic assignment by Frank-Wolfe algorithm.  This version uses MPI
 * to run a set of different inputs (changes to base net links file)
 * on multiple processors/nodes, each
 * of which runs the shortest path computations (d'Esopo-Pape
 * algorithm) in parallel using POSIX threads (one thread per source
 * node, up to max number of threads speicfieid to use).
 *
 *
 * There are many descriptions of this algorithm, a simple one for
 * example is in Sheffi, Y. 1985 "Urban Transportation Networks:
 * Equilibribum Analysis with Mathematical Programming Methods"
 * Prentice-Hall (out of print, available online from
 * http://web.mit.edu/sheffi/www/urbanTransportation.html
 *
 *   Usage: tap_frankwolfe_mpi [-n threads] [-w flowsfilename] [-c]
 *                                 [-i iter] [-r relgap]  [-p|-s|-l changelist]
 *                                 [-t tripmodfilenme]
 *                                 netfilename  demandfilename
 *                                 netmodfilename
 *                                 flowoutputbasefilename
 *
 *     netfilename    is name of the net file defining node and links
 *     demandfilename is name of Origin-Destination demand file
 *     netmodfilename  is name of net modification file defining changes
 *     flowoutputbasefilename is base name of flows output files
 *                        (network modification id is added for each entry
 *                         in netmodfilename for that modified network's
 *                         flow output)
 *     -n threads       : run each MPI process (each entry in netmodfilename)
 *                        with that many worker threads
 *     -w flowsfilename : warmstart from flows stored in flowsfilename
 *                        which is the output link flows file format
 *     -i iterations    : terminate after this many iterations
 *     -r relgap        : terminate after relative gap is <= relgap.
 *     -p               : run with each pair of modifications from the 
 *                        netmodfilename rather than each single modification
 *                        WARNING: this can be a large number of
 *                        modifications since there are n*(n-1)/2 pairs
 *                        of modifications for n modifications in the file.
 *     -s               : run with each subset (of size 3 or larger)
 *                        of modifications from the netmodfilename rather
 *                        than each single modification. 
 *                        WARNING: even more so than -p this can obviously
 *                        be very large since there are 2^n -n*(n-1)/2 -n -1
 *                        subsets of size 3 or larger. In fact this is
 *                        restricted to n = 8 modifictions at most.
 *    -l changelist     : instead of each modification separtely, run with just
 *                        the supplied list (comma delimited) of modifications,
 *                        all together, i.e. only run TAP is solved, so this 
 *                        option actually results in no parallelism (oher
 *                        than multithreading the single TAP solution).
 *    -t tripmodfilename: parse modifications to the O-D demand data form
 *                        tripmodfilename and run with the modified trip data
 *    -c                : continue from a previous interrupted run. This
 *                        is different from -w and does NOT imply it. This
 *                        option will check for each flow output name, and if 
 *                        already exists (or with .gz or .Z or .bz2 suffix),
 *                        will skip that modification, assuming it was already
 *                        run.
 *
 * WARNING: the flows output files constructed by appending the change
 * ids from the netmodfile to the flowoutputbasefilename will be 
 * overwritten if they exist.
 *
 * Note that becuase this uses both MPI and POSIX threads (pthreads),
 * it should be set up to run so that each MPI process is on a node
 * with enough cores to run the number of threads specified e.g.
 * run one MPI process on a single 8-core node with 8 threads or
 * two MPI processoes on a single 8-core node each with 4 threads, etc.
 * This is accomplished with the mpirun and PBS (where used, or similar)
 * job submission options, and the -n option on this program.
 * For example, to run 3 MPI processes (3 different inputs) each with
 * 4 threads, you would use PBS attributes such as:
 * #PBS -l nodes=3:ppn=4
 * and -np 3 -npersocket 1 on mpirun
 * (this is assuming a 'socket' is a quadcore CPU as for example on
 * tango.vpac.org)
 * to ensure there are 3 nodes each with 4 cores (of course this might
 * end up having two MPI slots  on a single  node, but that is OK
 * as there are 2 quad-core CPUs on each), and use -n 4 option
 * on tap_frankwolfe_mpi (since the default will be to use all "online"
 * cores which the sysconf() returns as 8 on an 8 core machine, regardles
 * of PBS attributes).
 *
 * MPI POSIX threads support (MPI_THREAD_MULTIPLE in OpenMPI) is not
 * required as we only do MPI in the main thread of each process, never
 * in other threads.
 *
 *   Example usage:
 *   
 * mpirun -np 3 -npersocket 1  ./tap_frankwolfe_mpi  -n 4              \
 * ~/traffic_assignment/trunk/testdata/SiouxFalls/SiouxFalls_net.txt   \
 * ~/traffic_assignment/trunk/testdata/SiouxFalls/SiouxFalls_trips.txt \
 * SiouxFalls_mods.txt SiouxFalls_flows
 *
 *   
 * See tap_frankwolfe_pthread.c for details of signal handling.
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
#include <sys/time.h>
#include <limits.h> /* PATH_MAX */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <mpi.h>

#include "utils.h"
#include "parsetapfiles.h"
#include "tap_functions.h"
#include "tap_frankwolfe_pthread.h"
   
/*****************************************************************************
 *
 * constants and typedefs
 *
 ****************************************************************************/

const int DEFAULT_ITERATION_LIMIT = 10000;

typedef struct int_pair_s
{
  int i;
  int j;
} int_pair_t;


#define MAXBIT 16 /* limit on number of mods for power set (-s) option */
static unsigned short bit[MAXBIT]; /* powers of 2 for powerset bit testijng */
#define set_bit(e,W) {W |= bit[e];}
#define test_bit(e,W) ((W) & bit[e])


   
/*****************************************************************************
 *
 * local functions
 *
 ****************************************************************************/

/*
 * count_bits_set() - count  number of bits set in a char
 *
 * Parameters:
 *     char to count bits set in
 *
 * Return value:
 *     number of bits set in the char
 *
 * http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetKernighan
 */
static unsigned char count_bits_set(unsigned short v)
{
  unsigned char c; /* accumulates total bits set in v */
  for (c = 0; v; c++)
    v &= v-1; /* clear LSB set */
  return c;
}

/* 
 * mod_entry_compar() - qsort comparison function for net_mod_t entries
 * 
 * Compares by change_id
 *
 */
static int mod_entry_compar(const void *ent1, const void *ent2)
{
  const net_mod_t *e1 = (const net_mod_t *)ent1;
  const net_mod_t *e2 = (const net_mod_t *)ent2;
  
  return strncmp(e1->change_id, e2->change_id, MAX_CHANGE_ID_LEN);
}

/*
 * apply_modification() -  Apply a single modification to the network data
 *
 * Parameters:
 *    net (in/out) - network data with links
 *    curmod - modificatino to make to the links in net data
 *
 * Return value:
 *    None.
 */
static void apply_modification(net_data_t *net, net_mod_t *curmod)
{
  int k;

  switch (curmod->modtype)
  {
    case MOD_TYPE_ADD_E:
      if (!(net->links = (link_data_t *)realloc(net->links, 
                                (net->num_links+1)*sizeof(link_data_t))))
      {
        fprintf(stderr, "realloc net->links failed\n");
            MPI_Finalize();exit(1);
      }
      /* TODO should check for duplicate links or bad node numbers etc. */
      net->links[net->num_links] = curmod->mod_link;
      ++net->num_links;
      fprintf(stderr, "[add %s] initnode %ld termnode %ld  capacity %f\n", curmod->change_id, net->links[net->num_links-1].init_node, net->links[net->num_links-1].term_node, net->links[net->num_links-1].capacity); /* XXX */
      break;
      
    case MOD_TYPE_CHANGE_E:
      /* TODO sort net data and do binary search instead */
      for (k = 0; k < net->num_links; k++)
        if (net->links[k].init_node == curmod->mod_link.init_node &&
            net->links[k].term_node == curmod->mod_link.term_node)
          break;
      if (k >= net->num_links)
      {
        fprintf(stderr, "ERROR: link not found for change %s\n", 
                curmod->change_id);
        MPI_Finalize(); exit(1);
      }
      fprintf(stderr, "[change %s] initnode %ld termnode %ld old capacity %f new capacity %f\n", curmod->change_id, net->links[k].init_node, net->links[k].term_node, net->links[k].capacity, curmod->mod_link.capacity); /* XXX */
      net->links[k] = curmod->mod_link;
      break;
      
    default:
      fprintf(stderr, "bad modifiction type %d\n", curmod->modtype);
      MPI_Finalize(); exit(1);
  }
}

/*****************************************************************************
 *
 * Main
 *
 ****************************************************************************/

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s [-c]  [-n numthreads] [-w warmstart_flows_filename] [-p|-s|-l changelist]  [-t tripmodsfilename] netfilename demandfilename netmodsfilename flowoutputbasefilename\n"
          "  -n numthreads : number of threads to use (default %d)\n"
          "  -w warmsetart_flows_filename : warmstart from previous solution\n"
          "  -i iterations : terminate after this many iterations (default %d)\n"
          "  -r relgap     : terminate when relative gap <= relgap\n"
          "  -p            : run all pairs of modifications in netmodsfilename\n"
          "  -s            : run all subsets  (size > 2) of modifications in netmodsfilename\n"
          "  -l changelist : run just the single TAP on network with list of change ids applied\n"
          "  -t tripmodsfilename : modifiy O-D demand data\n"
          "  -c            : continue from interrupted run,  do not redo existing flow output files\n"
          , progname, num_cores, DEFAULT_ITERATION_LIMIT );
  fprintf(stderr, "WARNING: the flows output files will be overwritten if they exist\n");
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
  int  numtasks, rank, rc; 
  double total_cost;
  int num_mods;
  int num_distinct_mods;
  char *modfile_name;
  FILE *modfile_fp;
  char *flow_output_basename;
  net_mod_t *mods;
  int i,k,j;
  int mod_i;
  char *prev_change_id;
  net_mod_t **change_id_index;
  net_mod_t *mod;
  net_mod_t *curmod;
  char myname[MPI_MAX_PROCESSOR_NAME]; int mynamelen;
  net_data_t net;
  struct timeval start_timeval, end_timeval, elapsed_timeval;
  int etime  = 0;
  char flow_output_filename[PATH_MAX+1];
  FILE *flow_output_fp;
  double relgap,objvalue;
  int iterations;
  int run_pairwise = 0;
  int run_subsets = 0;
  int_pair_t *pairindex = NULL;
  unsigned short *bitvecindex = NULL;
  unsigned short bitvec;
  int total_distinct_mods;
  int num_subsets_lt3;
  char changeid[MAX_CHANGE_ID_LEN*8+8]; /* changeid or changeid1_and_changeid2  or changeid1+changeid2+changeid3+...*/
  time_t now;
  char *trip_mods_filename = NULL;
  FILE *trip_mods_fp = NULL;
  int trip_mods_mode = 0;
  int num_trip_mods = 0;
  trip_mod_t *trip_mods = NULL;
  demand_data_t **demands = NULL;
  long num_zones = 0;
  int skip_existing_output_files = 0; /* 1 for -c option */
  char checkfilename[PATH_MAX+1];
  struct stat statbuf;
  int skip_this_one = 0;
  int run_with_changelist = 0;
  char *changelist = NULL;
  char **changeid_list = NULL;

  net.links = NULL;
  
  /* set up bit array for bit testing for -s option */
  c = 1;
  for (i = 0; i < MAXBIT; i++) {
    bit[i] = c;
    c *= 2;
  }
   
  num_cores = get_num_cores();
  num_threads = num_cores;
  
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    fprintf (stderr, "Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Get_processor_name(myname, &mynamelen);

  while ((c = getopt(argc, argv, "n:w:i:r:pst:cl:")) != -1)
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

      case 'p': /* do all pairs of modifications not indiviual modifications */
        run_pairwise = 1;
        break;

      case 's': /* do all subsets of size >=3 not indivdual modifications */
        run_subsets = 1;
        break;

      case 't': /* modifi origin-dest demand data not network */
        trip_mods_filename = optarg;
        trip_mods_mode = 1;
        break;

      case 'c': /* continue interrupted run; skip existing output files */
        skip_existing_output_files = 1;
        break;

      case 'l': /* run with supplied list of changes only (comma delimited) */
        run_with_changelist = 1;
        changelist = optarg;
        break;

      default:
        usage(argv[0]);
        break;
    }
  }
  
  if (argc - optind != 4)
    usage(argv[0]);


  if (run_pairwise + run_subsets + run_with_changelist > 1)
  {
    fprintf(stderr, "-p and -s and -l changelist are all mutually exclusive\n");
    usage(argv[0]);
  }

  net_filename = argv[optind];
  demand_filename = argv[optind+1];
  modfile_name = argv[optind+2];
  flow_output_basename = argv[optind+3];

  if (target_relgap <= 0 && target_iterations <= 0)
    target_iterations = DEFAULT_ITERATION_LIMIT;


  if (trip_mods_mode)
  {
    if (!(trip_mods_fp = fopen(trip_mods_filename, "r")))
    {
      fprintf(stderr, "error opening trips modifications file %s: %s\n",
              trip_mods_filename, strerror(errno));
      MPI_Finalize(); exit(1);
    }
    if (parse_trip_mod_file(trip_mods_fp, &trip_mods, &num_trip_mods) != 0)
    {
      fprintf(stderr, "parse trip mods file %s failed\n", trip_mods_filename);
      MPI_Finalize(); exit(1);
    }
    fclose(trip_mods_fp);
  }

  /*
   * Each MPI process parses the network modifications file and
   * selects for itself a fair share of the modifications to process,
   * writing the flows at the end for each modification to a file
   * with the changeId for that modification in the filename.
   */
  if (!(modfile_fp = fopen(modfile_name, "r")))
  {
    fprintf(stderr, "error opening modifications file %s: %s\n",
            modfile_name, strerror(errno));
    MPI_Finalize(); exit(1);
  }
  
  if (parse_net_mod_file(modfile_fp, &mods, &num_mods) != 0)
  {
    fprintf(stderr, "parse net mods file %s failed\n", modfile_name);
    MPI_Finalize(); exit(1);
  }
  fclose(modfile_fp);

  /* find the number of distinct change ids (a single change id
     can be composed of multiple changes that are considred indivisible).
     Also build quick lookup list of pointers to first entry for each
     distinct change id */
  if (!(change_id_index = (net_mod_t **)malloc(num_mods * sizeof(net_mod_t))))
  {
    fprintf(stderr, "malloc change_id_index failed\n");
    MPI_Finalize(); exit(1);
  }
  qsort(mods, num_mods, sizeof(net_mod_t), mod_entry_compar);
  num_distinct_mods = 0;
  prev_change_id = "";
  for (i = 0; i < num_mods; i++)
  {
    if (strcmp(mods[i].change_id, prev_change_id) != 0)
    {
      change_id_index[num_distinct_mods++] = &mods[i];
      prev_change_id = mods[i].change_id;
    }
  }

/*  fprintf(stderr, "%d distinct mods\n",num_distinct_mods); */

  if (run_pairwise)
  {
    /* build index mapping sequential number to a (i,j) pair
       ofr change indices so we can just do a single loop below.
       TODO
       There is probably a better way to do this than building this lookup
       table, e.g.
       http://www.cs.berkeley.edu/~wkahan/Math55/pairs.pdf 
       but this is easy and it works - we don't care so much about
       the waste of memory since for large values of num_distinct_mods
       the number of things we have to process becomes impractical anyway
    */
    total_distinct_mods =  num_distinct_mods*(num_distinct_mods-1)/2;
    if (rank == 0)
      fprintf(stderr, "running pairwise modifications: %d total\n", 
              total_distinct_mods);
    if (!(pairindex = malloc(total_distinct_mods*sizeof(int_pair_t))))
    {
      fprintf(stderr, "malloc pairindex failed\n");
      MPI_Finalize(); exit(1);
    }
    k = 0;
    for (i = 0; i < num_distinct_mods; i++)
    {
      for (j = i+1; j < num_distinct_mods; j++)
      {
        pairindex[k].i = i;
        pairindex[k].j = j;
        ++k;
      }
    }
  }
  else if (run_subsets)
  {
    /* We can just use a bit vector to select each subset, by counting
       an 8-bit integer from 0 to n and each bit indicating whether
       the correspdoning changeid is in the subset or not.
       But we also want to just skip over any subset of size less than 3
       (the empty set, single changes, pairwise changes) as these are 
       either trivial or done with other options (default, -p).
    */
    if (num_distinct_mods > MAXBIT)
    {
      fprintf(stderr, "error: maximum of %d modifications allowed for -s "
              "option but %d in file %s\n", MAXBIT, num_distinct_mods,
              modfile_name);
      MPI_Finalize(); exit(1);
    }
    num_subsets_lt3 = num_distinct_mods*(num_distinct_mods-1)/2 /* size 2 */
      + num_distinct_mods /*size 1*/ + 1 /*size 0*/;
    total_distinct_mods = bit[num_distinct_mods-1]*2 /* 2^num_distinct_mods */
      - num_subsets_lt3; /* number of subsets of size less than 3 */
    if (!(bitvecindex = malloc(total_distinct_mods*sizeof(unsigned short))))
    {
      fprintf(stderr, "malloc bitvecindex failed\n");
      MPI_Finalize(); exit(1);
    }
    if (rank == 0)
      fprintf(stderr, "running power set of modifications: %d total\n",
              total_distinct_mods);
    /* build index mapping sequential number to bitvector for subset,
       excluding those subsets with size < 3 */
    bitvec = 0x7; /* 7 is first with more than 2 bits set */
    k = 0;
    while (k < total_distinct_mods)
    {
      if (count_bits_set(bitvec) >= 3)
        bitvecindex[k++] = bitvec;
      ++bitvec;
    }
  }
  else if (run_with_changelist)
  {
    char *saveptr = NULL;
    char *token;
    int l;

    total_distinct_mods = 1; /* the list of changes is made into one combined change */
    if (!(changeid_list = (char **)calloc(num_distinct_mods , sizeof(char*))))
    {
      fprintf(stderr, "malloc changeid_list failed\n");
      MPI_Finalize(); exit(1);
    }
    token = strtok_r(changelist, ",", &saveptr);
    changeid_list[0] = token;
    for (l = 1; token != NULL && l < num_distinct_mods; l++)
    {
      token = strtok_r(NULL, ",", &saveptr);
      changeid_list[l] = token;
    }
    if (rank == 0)
    {
      fprintf(stderr, "MPI rank %d (%s) processing change: ", rank, myname);
      l = 0;
      while (changeid_list[l])
      {
        fputs(changeid_list[l], stderr);
        if (changeid_list[l+1])
          fputs(" + ", stderr);
        ++l;
      }
      fputs("\n", stderr);
    }
    else /* shoudl only have 1 MPI process for -l mode, no parallelism here */
      fprintf(stderr, "MPI rank %d (%s) idle for -l mode\n", rank, myname);
      
  }
  else
    total_distinct_mods = num_distinct_mods;

  /* each MPI process will process ceil(total_distinct_mods / numtasks)
     distinct change_ids */
  for (mod_i = rank; mod_i < total_distinct_mods; mod_i += numtasks)
  {
    
    if (!(net_fp = fopen(net_filename, "r")))
    {
      fprintf(stderr, "error opening net file %s: %s\n", 
              net_filename, strerror(errno));
      MPI_Finalize(); exit(1);
    }
    
    if (!(demand_fp = fopen(demand_filename, "r")))
    {
      fprintf(stderr, "error opening trips file %s: %s\n",
              demand_filename, strerror(errno));
      MPI_Finalize(); exit(1);
    }
    
    if (warm_start_mode)
    {
      /* we ahve the 'warm start' option to parse link volumes (and costs)
         from a flows file (output of this program also) */
      if (!(flows_input_fp = fopen(flows_input_filename, "r")))
      {
        fprintf(stderr, "error opening flows input file %s: %s\n",
                flows_input_filename, strerror(errno));
        MPI_Finalize(); exit(1);
      }
    }

    if (parse_net_file(net_fp, &net) != 0)
    {
      fprintf(stderr, "error parsing net file\n");
      MPI_Finalize(); exit(1);
    }
    fclose(net_fp);
    

    if (run_pairwise)
    {
      strncpy(changeid, change_id_index[pairindex[mod_i].i]->change_id, sizeof(changeid)-1);
      strncat(changeid, "_and_", sizeof(changeid)-1);
      strncat(changeid, change_id_index[pairindex[mod_i].j]->change_id, sizeof(changeid)-1);
    }
    else if (run_subsets)
    {
      changeid[0] = '\0';
      for (k = 0; k < num_distinct_mods; k++)
      {
        if (test_bit(k, bitvecindex[mod_i]))
        {
          strncat(changeid, change_id_index[k]->change_id, sizeof(changeid)-1);
          strncat(changeid, "+", sizeof(changeid)-1);
        }
      }
      changeid[strlen(changeid)-1] = '\0'; /* remove '+' on end */
    }
    else if (run_with_changelist)
    {
      changeid[0] = '\0';
      k = 0;
      while (changeid_list[k])
      {
        strncat(changeid, changeid_list[k], sizeof(changeid)-1);
        strncat(changeid, "+", sizeof(changeid)-1);
        ++k;
      }
      changeid[strlen(changeid)-1] = '\0'; /* remove '+' on end */
    }
    else
      strncpy(changeid, change_id_index[mod_i]->change_id, MAX_CHANGE_ID_LEN);

    fprintf(stderr, "MPI rank %d (%s) processing id %s\n", rank,
            myname, changeid);
    

    if (run_pairwise)
    {
      /* make the changes to the network from the current pair of change ids
         (there may be multiple single changes for the same change id) */
      for (curmod = mod = change_id_index[pairindex[mod_i].i]; 
           strncmp(curmod->change_id, mod->change_id, MAX_CHANGE_ID_LEN) == 0;
           curmod++)
        apply_modification(&net, curmod);
      for (curmod = mod = change_id_index[pairindex[mod_i].j]; 
           strncmp(curmod->change_id, mod->change_id, MAX_CHANGE_ID_LEN) == 0;
           curmod++)
        apply_modification(&net, curmod);
    }
    else if (run_subsets)
    {
      for (k = 0; k < num_distinct_mods; k++)
      {
        if (test_bit(k, bitvecindex[mod_i]))
        {
          for (curmod = mod = change_id_index[k]; 
               strncmp(curmod->change_id,
                       mod->change_id, MAX_CHANGE_ID_LEN) == 0;
               curmod++)
          {
            apply_modification(&net, curmod);
          }
        }
      }
    }
    else if (run_with_changelist)
    {
      /* make all the changes in the list to th enetwork */
      /* (there may be multiple single changes for the same change id) */
      /* TODO this is quadratic, should build an index to make it linear */
      /* TODO should also check for duplicate mods in list, etc. */
      k = 0;
      while (changeid_list[k])
      {
        int foundit = 0 ;
        for (i = 0; i < num_distinct_mods; i++)
        {
          if (strncmp(change_id_index[i], changeid_list[k], MAX_CHANGE_ID_LEN) == 0)
          {
            foundit = 1;
            for(curmod = mod = change_id_index[i];
                strncmp(curmod->change_id, mod->change_id, MAX_CHANGE_ID_LEN) == 0;
                curmod++)
            {
              apply_modification(&net ,curmod);
            }
          }
        }
        if (!foundit)
        {
          fprintf(stderr, "ERROR: change %s not found\n", changeid_list[k]);
          MPI_Finalize(); exit(1);
        }
        ++k;
      }
    }
    else
    {
      /* make the changes to the network from the current change id
         (there may be multiple single changes for the same change id) */
      for (curmod = mod = change_id_index[mod_i]; 
           strncmp(curmod->change_id, mod->change_id, MAX_CHANGE_ID_LEN) == 0;
           curmod++)
      {
        apply_modification(&net, curmod);
      }
    }

    /* sort link data so matches saved flow data if warmstart used */
    /* TODO can remove this if already sorted and insertion sort used
       for new modificatino links if TODO above is implemnted */
    qsort(net.links, net.num_links, sizeof(link_data_t), link_data_compar);

    strncpy(flow_output_filename, flow_output_basename, PATH_MAX);
    strncat(flow_output_filename, changeid, PATH_MAX);
    strncat(flow_output_filename, ".txt", PATH_MAX);
    if (skip_existing_output_files)
    {
      strncpy(checkfilename, flow_output_filename, PATH_MAX);
      skip_this_one = 0;
      if (stat(flow_output_filename, &statbuf) == 0)
        skip_this_one = 1;
      else 
      {
          strncpy(checkfilename, flow_output_filename, PATH_MAX);
          strncat(checkfilename,  ".gz", PATH_MAX);
          if (stat(checkfilename, &statbuf) == 0)
              skip_this_one = 1;
          else 
          {
              strncpy(checkfilename, flow_output_filename, PATH_MAX);
              strncat(checkfilename,  ".Z", PATH_MAX);
              if (stat(checkfilename, &statbuf) == 0)
                  skip_this_one = 1;
              else 
              {
                  strncpy(checkfilename, flow_output_filename, PATH_MAX);
                  strncat(checkfilename,  ".bz2", PATH_MAX);
                  if (stat(checkfilename, &statbuf) == 0)
                      skip_this_one = 1;
              }
          }
       }
       if (skip_this_one)
       {
           fprintf(stderr, "MPI rank %d (%s) flow output file %s exists, skipping\n",
                   rank, myname,  checkfilename);
           continue; /* skip this modificatin, do next one */
       }
    }

    if (!(flow_output_fp = fopen(flow_output_filename, "w")))
    {
      fprintf(stderr, "ERROR: cannot open %s for writing (%s)\n",
              flow_output_filename, strerror(errno));
      fprintf(stderr, "    Flow output for change id %s will be written to stdout\n", changeid);
      flow_output_fp = stdout;
    }

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

    gettimeofday(&start_timeval, NULL);
    if (tap_frankwolfe(&net, demands,
                       warm_start_mode, 
                       flows_input_fp,
                       target_iterations, target_relgap, flow_output_fp,
                       &total_cost, &objvalue, &iterations, &relgap, 0) != 0)
    {
      fprintf(stderr, "tap_frankwolfe() failed\n");
      MPI_Finalize(); exit(1);
    }
    gettimeofday(&end_timeval, NULL);
    timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
    etime = 1000*elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
    fprintf(stderr, "MPI rank %d (%s) id %s iterations %d objective value %.15f relative gap %.15f time %.2fs\n", rank,
            myname, changeid, iterations, objvalue, relgap, (double)etime/1000);

    fprintf(stderr, "id = %s VHT = %f\n",  changeid, total_cost);
    fprintf(flow_output_fp, "~ iterations = %d objvalue = %.15f relgap = %.15f\n", iterations, objvalue, relgap);
    fprintf(flow_output_fp, "~ totalcost = %.15f\n", total_cost);
    fprintf(flow_output_fp, "~ mpirank = %d machine = %s threads = %d time = %.2fs\n", rank, myname, num_threads,  (double)etime/1000);
    fprintf(flow_output_fp, "~ generated by: ");
    for (i = 0; i < argc; i++)
      fprintf(flow_output_fp, "%s ", argv[i]);
    fprintf(flow_output_fp, "\n");
    now = time(NULL);
    fprintf(flow_output_fp, "~ on: %s", ctime(&now));
    fprintf(flow_output_fp, "~ changeid = %s\n", changeid);

    fclose(flow_output_fp);
  }


  /* TODO use MPI to send results to rank 0 process to compute
     correlations between pairs of modifications etc.
  */

  free(net.links);
  free(change_id_index);
  free(pairindex); 
  free(bitvecindex);
  free(changeid_list);
  MPI_Finalize();
  exit(0);
}

