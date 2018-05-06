/*****************************************************************************
 * 
 * File:    sssp.cu
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: sssp.cu 855 2011-11-15 02:47:59Z astivala $
 *
 * Single-source shortest path implementation using CUDA, double precision
 * 
 ****************************************************************************/

/* FIXME -  should either fix the start_nodes vector code to actually allow arbitrary start nodes or remove it and hardcode to have start nodes 1..num_zones  which it currently assumes (with hack to add 1 s for not 0.. num_zones-1) */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <cutil_inline.h>      /* CUDA SDK */

#include "sssp.h"
#include "sssp_gold.h"
#include "sssp_pape.h"
#include "pape_cuda_host.h"
#include "bellmanford_spblas.h"
#include "bellmanford_cusp.h"

#undef DUMP_RESULTS
#undef TRACEBACK
#undef ONLY_RUN_CUDA
#define RUN_SLOW_BELLMANFORD
#define RUN_PAPE_LLL_SAME_MEMORY
#undef RUN_SPARSE_BLAS_BELLMAN_FORD /* FIXME there seems to be a memory leak in SparseBLAS library, don't use this */
#undef RUN_CUSP_BELLMAN_FORD /* this is too slow to be practical */
#undef DEBUG

/*****************************************************************************
 *
 * constants
 *
 ****************************************************************************/

const double EPS = 1e-08;  // error toelrance

/*****************************************************************************
 *
 * globals
 *
 ****************************************************************************/
bool using_fermi_architecture = false; // using a Fermi architecture GPU

/*****************************************************************************
 *
 * Local functions
 *
 ****************************************************************************/

//Align a to nearest higher multiple of b
extern "C" long iAlignUp(int a, long b){
    return ((a % b) != 0) ?  (a - a % b + b) : a;
}

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
void dump_packed_arrays(const long Va[], const long Ea[], const double Wa[],
                        long num_nodes, long num_edges)
{
  long i;
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


/*
 * verify_results() -compare two doubleing point arrays 
 *
 * Parameters:
 *    test -  results to check
 *    ref - reference versino to compare against 
 *    n - length of arrays to compare
 *
 * Return value:
 *    0 if OK else -1
 */
static long verify_results(double *test, double *ref, long n)
{
  long rc = 0;
  for (long i = 0; i < n; i++)
    if (fabs(ref[i] - test[i]) > EPS)
    {
      fprintf(stdout, "i = %ld ref = %f test = %f\n", i, ref[i], test[i]);
      rc = -1;
    }
  return rc;
}



/* Retreive the mincost path from i to j implied by the
   distance and predecessor matrices
   NB this path is stored backwards (from j to i) */
int Min_getpath(double *distances, long *predecessors, long num_nodes,
                long num_start_nodes, long i, long j, long **path, 
                long *pathlen, bool sourcemajor)
{
  long idx = sourcemajor ? i * num_nodes + j : j * num_start_nodes + i;
  if (distances[idx] == FLOATINF)
    return 0; // no path
  long k = j;
  while (k != i)
  {
    assert(*pathlen < num_nodes);
    **path = k;
    (*path)++;
    (*pathlen)++;
    long k_idx = sourcemajor ? i * num_nodes + k : k * num_start_nodes + i;
    k = predecessors[k_idx];
  }
  **path = i;
  (*path)++;
  (*pathlen)++;
  return 0;
}

/*
 * verify_path() - verify that a traceback path and its cost is valid
 *
 * Parameters:
 *    path - stored in reverse from destinatino to source
 *    pathlen - length of path (number of nodes inc. dest and source)
 *    cost - purported cost of path
 *    Va[num_nodes] -array of indices to head of each adj list in Ea 
 *    Ea[num_edges] -each entry is 'to' node in list for 'from' node
 *    Wa[num_edges] -each entry is cost of corresponding Ea entry
 *    num_nodes  -number of nodes
 *    num_edges  -number of edges
 *
 * Return value:
 *   0 if ok else -1
 */
static long verify_path(long path[], long pathlen, double cost,
                       long Va[], long Ea[], double Wa[], long num_nodes,
                       long num_edges)
{
  if (cost == FLOATINF && pathlen < 1)
    return 0; //no path
  double acost = 0;
  for (long k = pathlen-1; k > 0; k--)
  {
    if (path[k] > num_nodes)
    {
      fprintf(stderr, "bad path: %ld > num_nodes at %ld\n", path[k], k);
      return -1;
    }
    if (path[k-1] > num_nodes)
    {
      fprintf(stderr, "bad path: %ld > num_nodes at %ld\n", path[k+1], k+1);
      return -1;
    }

    long start = Va[path[k]];
    long end = Va[path[k]+1];
    long i;
    for (i = start; i < end && Ea[i] != path[k-1]; i++)
      /*nothing*/;
    if (Ea[i] != path[k-1])
    {
      fprintf(stderr, "bad path: no edge from %ld to %ld\n", path[k],
              path[k-1]);
      return -1;
    }
    acost += Wa[i];
  }
  if (fabs(acost - cost) > EPS)
  {
    fprintf(stderr, "bad path: cost %f != %f\n", acost, cost);
    return -1;
  }
  return 0;
}


/*
 * check_paths() - verify that traceback paths and their costs are valid
 *
 * Parameters:
 *    Va[num_nodes] -array of indices to head of each adj list in Ea 
 *    Ea[num_edges] -each entry is 'to' node in list for 'from' node
 *    Wa[num_edges] -each entry is cost of corresponding Ea entry
 *    num_nodes  -number of nodes
 *    num_edges  -number of edges
 *    num_start_nodes - number of soruce nodes
 *    distances  - shortest path distance matrix
 *    predecessor - predcessor matrix
 *    sourcemajor - true if rows are each source, else rows each dest vertex
 *
 * Return value:
 *   0 if ok else -1
 */
static long check_paths(long Va[], long Ea[], double Wa[], long num_nodes,
                       long num_edges, long num_start_nodes,
                       double distances[],
                       long predecessors[], bool sourcemajor)
{

  /* check paths */
  bool pathtest_failed = false;
  long *thepath  = new long[num_nodes];
  for (long i = 0; i < num_start_nodes; i++)
  {
    for (long j = 0; j < num_nodes; j++)
    {
      long *path = thepath;
      long **pathptr = &path;
      long pathlen = 0;
      Min_getpath(distances, predecessors, num_nodes, num_start_nodes,
                  i, j, pathptr, &pathlen, sourcemajor);
      long idx = sourcemajor ? i * num_nodes + j : j * num_start_nodes + i;
      if (verify_path(thepath, pathlen, 
                      distances[idx],
                      Va, Ea, Wa, num_nodes,
                      num_edges) != 0)
      {
        printf("verify path FAILED source %ld dest %ld\n", i, j);
        pathtest_failed = true;
      }
    }
  }
  if (!pathtest_failed)
    printf("verify %ld paths ok\n", num_start_nodes*num_nodes);
  delete[] thepath;
  if (pathtest_failed)
    return -1;
  else
    return 0;
}


/*
 * read_adjlist() - read adjency list from file
 *
 * Parameters:
 *    fp - open (read) filepointer to read from
 *    adjlist (OUT) - newly allocated adjacenty list 
 *    num_nodes (OUT) - numer of nodes (max node number)
 *
 * Return value:
 *    number of entries read (length of adjlist); -1 on error
 *
 * Reads the adjcancy list from file in format used by the apsp program
 * ie.:
 *
 *    tab-delimited lines
 *    first row is
 *    x  y  z
 *    where x,y are dimensions (x=y) and nz is number of nonzero entries
 *    each subsequent row is
 *
 *    i  j  c
 *
 *    where  i and j are (1-based) node nubmers and c is cost from i to j
 *
 *   e.g.
 *   8       8       10
 *   11      12      0.643311
 *   12      13      0.775850
 *   13      14      0.000036
 *   14      15      0.079787
 *   15      16      0.483686
 *   16      17      0.883406
 *   17      18      0.281150
 *   18      13      0.055000
 *   12      14      0.093787
 *
 * Note x and y not actually used here
 *
 * NB tht input is 1-baed but we subtract 1 to be zero-baed internally
 * 
 */  
static long read_adjlist(FILE *fp, adjlist_entry_t **adjlist, long *num_nodes)
{
  long x,y,nnz,n;
  adjlist_entry_t *entry;
  long max_node_num = 0;

  if (fscanf(fp, "%d\t%d\t%d\n", &x, &y, &nnz) != 3)
  {
    fprintf(stderr, "bad first line in input file\n");
    return -1;
  }
  if (!(*adjlist = (adjlist_entry_t *)malloc(nnz*sizeof(adjlist_entry_t))))
  {
    fprintf(stderr, "malloc adjlist failed\n");
    return -1;
  }
  n = 0;
  entry = *adjlist;
  while (!feof(fp))
  {
    if (fscanf(fp, "%d\t%d\t%lf\n",&entry->from, &entry->to, &entry->cost)!= 3) /* NB scanf() REQUIRES %lf for double to work, while in printf() %f works for double */
    {
      fprintf(stderr, "error in input file at %dth entry\n", n);
      return -1;
    }
    if (++n > nnz)
    {
      fprintf(stderr, "too many entries in input file (cutoff at %d)\n",n);
      return -1;
    }
    if (entry->from > max_node_num || entry->to > max_node_num)
      max_node_num = MAX(entry->from, entry->to);

    // convert to 0-based node numbers
    --entry->from;
    --entry->to;

    ++entry;
  }
  *num_nodes = max_node_num;
  return n;
}

/* 
 * adjlist_entry_compar() - qsort comparison function for adlist entries
 * 
 * Compares by 'from' node number first  then by 'to' node number if equal
 *
 */
static int adjlist_entry_compar(const void *ent1, const void *ent2)
{
  const adjlist_entry_t *e1 = (const adjlist_entry_t *)ent1;
  const adjlist_entry_t *e2 = (const adjlist_entry_t *)ent2;
  
  if (e1->from < e2->from)
    return -1;
  else if(e1->from > e2->from)
    return 1;
  else
    return e1->to < e2->to ? -1 : (e1->to > e2->to ? 1 : 0);
}

/*
 * adjlist_to_packed_arrays() - convert adjlist struct array to packed arrays
 *
 * Packed array format is illustrated in the referenced papers. It is like
 * like the packed column ("Harwell-Boeing") format used for sparse
 * matrices in some FORTRAN linear algebra routines.
 * For each node i, Va[i] is the first and Va[i+1]-1 the last index into
 * Ea containing the adjacent nodes to i and Wa giving the costs of
 * those edges from node i.
 *
 * Parameters:
 *        adjlist (In/Out) - array of adjlist entries (from,to,cost)
 *        num_nodes - number of nodes
 *        num_edges -  number of edges (length of adjlist)
 *        Va[num_nodes+1] (OUT)-array of indices to head of each adj list in Ea 
 *        Ea[num_edges] (OUT)-each entry is 'to' node in list for 'from' node
 *        Wa[num_edges] (OUT)-each entry is cost of corresponding Ea entry
 *
 * Return value:
 *     None.
 *
 * The Va, Ea and Wa arrays must be already allocated to the correct sizes
 * by the caller (num_nodes, num_edges and num_edges respectively).
 * The adjlist input is sorted in-place.
 */
static void adjlist_to_packed_arrays(adjlist_entry_t adjlist[],
                                     long num_nodes,
                                     long num_edges,
                                     long Va[],
                                     long Ea[],
                                     double Wa[])
{
  // sort by 'from' node ascending and within that by 'to' node ascending
  qsort(adjlist, num_edges, sizeof(adjlist_entry_t), adjlist_entry_compar);
  long v = 0;  // index into Va
  long e = 0;  // index into Ea and Wa and adjlist

  for (v = 0; v < num_nodes; v++)
  {
    Va[v] = e;
    while (e < num_edges && adjlist[e].from == v)
    {
      Ea[e] = adjlist[e].to;
      Wa[e] = adjlist[e].cost;
      e++;
    }
  }
  Va[num_nodes] = e;
}

/*
 * choose_gpu() - choose the GPU to use
 *
 * Parameters: 
 *      major_req - required (minimum) major compute capability
 *      minor_req - required (minimum) minor compute capabililty
 * Return value: 0 if OK, nonzero if no suitable device found
 *
 * Sets the global using_fermi_architecture and calls cudaSetDevice()
 * to set the GPU to use.
 */
static int choose_gpu(int major_req, int minor_req)
{
/*
    long devnum = cutGetMaxGflopsDeviceId();
    fprintf(stderr, "using max gflops device %d: ", devnum);
*/
    /* If there is a compute capability 2 device ("Fermi"
       architecture) (or higher) then use that
    */

    int devnum, deviceCount, gflops,max_gflops=0, sel_devnum=-1;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
    {
      fprintf(stderr, "There is no device supporting CUDA.\n");
      return 1;
    }
    for (devnum = 0; devnum < deviceCount; devnum++)
    {  
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, devnum);
      if (devnum == 0 && 
          deviceProp.major == 9999 && deviceProp.minor == 9999)
      {
        fprintf(stderr, "There is no device supporting CUDA.\n");
        return 1;
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
        using_fermi_architecture = true;
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
          using_fermi_architecture = false;
        }
      }
    }
    if (sel_devnum < 0)
    {
      fprintf(stderr,
              "there are no CUDA devices of required compute capability\n");
      return 2;
    }

    fprintf(stdout, "using device %d: ", sel_devnum);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, sel_devnum);
    fprintf(stdout, "%s\n", deviceProp.name);
    cudaSetDevice( sel_devnum );
    return 0;
}

#ifdef TRACEBACK
/*
 * printpaths() - print path given cost and predecesor vectors
 *
 * Parameters:
 *      fp - open (write) filepointer to write to
 *      cost - cost vectors from source to each node
 *      predecessor - predecessor node vectors from source to each node
 *      source  - source node
 *      num_nodes - number of nodes 
 *      num_start_nodes -number of source nodes
 *      sourcemajor - true if rows are each source, else rows each dest vertex
 *
 * Return value:
 *      None.
 *
 */
static void printpaths(FILE *fp, double cost[], long predecessor[], long source,
                       long num_nodes, long num_start_nodes, bool sourcemajor)
{
  for (long dest = 0; dest < num_nodes; dest++)
  {
    long idx = sourcemajor ? source*num_nodes+dest : dest*num_start_nodes+source;
      
    if (cost[idx] == FLOATINF)
    {
      fprintf(fp, "%ld -> %ld: no path\n", source+1, dest+1);
      continue;
    }
    fprintf(fp, "%ld -> %ld (%f): ", source+1, dest+1, cost[idx]);
    long k = dest;
    long k_idx = sourcemajor ? source*num_nodes+k : k*num_start_nodes+source;
    while (k != source)
    {
      fprintf(fp, "%ld <- ", k+1);
      k = predecessor[k_idx];
      k_idx = sourcemajor ? source*num_nodes+k : k*num_start_nodes+source;
    }
    fprintf(fp, "%ld \n", source+1);
  }
}
#endif /* TRACEBACK */


/*****************************************************************************
 *
 * Main
 *
 ****************************************************************************/

void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s inputfile\n", progname);
  exit(1);
}

int main(int argc, char *argv[])
{
  unsigned int hTimer;
  double runtime;
  double *distances;
  long *predecessors;
  bool failed = false, anyfailed = false;
  const long NUM_RUNS = 10; // run each test this many times to avg timings
  long *work_queue_next;

  if (argc != 2)
    usage(argv[0]);

  /* Read adjacney list from input file */

  char *inputfilename = argv[1];
  FILE *fp = fopen(inputfilename, "r");
  if (!fp)
  {
    fprintf(stderr, "open %s (read) failed\n", inputfilename);
    exit(1);
  }
  adjlist_entry_t *adjlist; long num_nodes = 0;
  long  num_edges = read_adjlist(fp, &adjlist, &num_nodes);
  if (num_edges < 0)
  {
    fprintf(stderr, "reading input file %s failed\n", inputfilename);
    exit(1);
  }

  /* Initialize GPU */
  bool use_cuda = 1;
  // need at least compute capability 1.3 for double precision f.p.
  if (choose_gpu(1,3) != 0) /* calls cudaSetDevice() */
  {
    fprintf(stderr, "WARNING: CUDA implementation disabled as no suitable device\n");
    use_cuda = 0;
  }


  /* first we have to convert adjacnecy lists to packed array format */

  // pad the arrays out to a multiple of BLOCK_SIZE
//  long alloc_num_nodes =  iAlignUp(num_nodes, BLOCK_SIZE);
  long alloc_num_nodes = num_nodes+1; //no longer need the above
  long *Va = new long[alloc_num_nodes];
  long *Ea = new long[num_edges];
  double *Wa = new double[num_edges];

  adjlist_to_packed_arrays(adjlist, num_nodes, num_edges, Va, Ea, Wa);

#ifdef DEBUG
  dump_packed_arrays(Va, Ea, Wa, num_nodes, num_edges);
#endif  

//  const long NUM_START_NODES =  147; // # of zones in winnipeg data // num_nodes; // / 2; //; // XXX
 const long NUM_START_NODES =  1790; // XXX # of zones in ChicagoRegional
//const long NUM_START_NODES =  2253; // XXX # of zones in Melbourne
//const long NUM_START_NODES = num_nodes / 2;
//const long NUM_START_NODES = 1; //XXX
  long num_start_nodes = NUM_START_NODES;

  assert(num_start_nodes <= num_nodes);

  long *start_nodes = new long[num_start_nodes];
  for (long i = 0; i < num_start_nodes; i++)
    start_nodes[i]  = i;

  double *gold_distances = new double[num_nodes*num_start_nodes];
  long *gold_predecessors = new long[num_nodes*num_start_nodes];


  fprintf(stdout, "input graph is %s with %ld nodes and %ld edges (%d KB in packed format)\n",
          inputfilename, num_nodes, num_edges,
          (num_nodes*sizeof(long)+num_edges*sizeof(long)+num_edges*sizeof(double))/1024);
  fprintf(stderr, "number of start nodes: %d\n", num_start_nodes);
  fprintf(stderr, "number of runs of each algorithm: %d\n", NUM_RUNS);
  int run_num = 0;
  double total_time = 0.0;
  double total_internal_time = 0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Runing Dijkstra algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    double internal_time = 0;
    /* Run on host using Boost Graph Library for verification */
    for (long i = 0; i < num_start_nodes; i++)
      internal_time += sssp_gold(adjlist, num_nodes, num_edges, start_nodes[i],
                gold_distances + i*num_nodes, gold_predecessors + i*num_nodes);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU Dijkstra algorithm internal  time %f ms\n",
             internal_time);
    fprintf(stdout,  "CPU Dijkstra algorithm execution time %f ms\n", runtime);
    total_time +=  runtime;
    total_internal_time += internal_time;
  }
  fprintf(stdout,  "CPU Dijkstra\talgorithm average time %f ms\n", 
          total_time / NUM_RUNS);
  fprintf(stdout,  "CPU Dijkstra\talgorithm average internal time %f ms\n", 
          total_internal_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("reference.paths", "w");
  for (long s = 0; s < num_start_nodes; s++)
    printpaths(fp, gold_distances, gold_predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */


  distances = new double[num_nodes*num_start_nodes];
  predecessors = new long[num_nodes*num_start_nodes];

#ifndef ONLY_RUN_CUDA
  long *queue_next;
  if (!(queue_next = (long *)malloc(num_nodes *sizeof(long))))
  {
    fprintf(stderr, "malloc queue failed\n");
    exit(1);
  }

  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Runing d'Esopo-Pape algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );


    for (long i = 0; i < num_start_nodes; i++)
      sssp_pape(Va, Ea, Wa, num_nodes, num_edges, start_nodes[i], 0,
                distances + i*num_nodes, predecessors + i*num_nodes,
                queue_next);
    /* reusing same distances/prdecessors vector rather than storing all results
       makes some difference, (eg 26s not 31s) but not the order of magnitude
       found in actual use ??? */
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU d'Esopo-Pape algorithm execution time %f ms\n", runtime);
    total_time +=  runtime;
    failed = false;
    if (run_num == 0)
    {
      /* compare results */
      for (long i = 0; i < num_start_nodes; i++)
      {
        if (verify_results(distances + i*num_nodes,
                           gold_distances + i*num_nodes, num_nodes) != 0)
        {
          fprintf(stdout, "test FAILED (i = %d)\n", i);
          failed = true; anyfailed = true;
        }
      }
      if (!failed)
        fprintf(stdout, "test %d costs ok\n", num_start_nodes *num_nodes);
      
      /* check paths */
      check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                  distances, predecessors, true);
    }
  }
  fprintf(stdout,  "CPU d'Esopo-Pape\talgorithm average time %f ms\n", 
          total_time / NUM_RUNS);

  free(queue_next); queue_next = NULL;

#ifdef TRACEBACK
  fp = fopen("test_pape.paths", "w");
  for (long s = 0; s < num_start_nodes; s++)
    printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */

  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Running d'Esopo-Pape LLL algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    for (long i = 0; i < num_start_nodes; i++)
    {
      if (!(work_queue_next = (long *)malloc(num_nodes *sizeof(long))))
      {
        fprintf(stderr, "malloc queue failed\n");
        exit(1);
      }

      sssp_pape_lll(Va, Ea, Wa, num_nodes, num_edges, start_nodes[i], 0,
                    distances + i*num_nodes, predecessors + i*num_nodes, work_queue_next);

      free(work_queue_next);
    }

    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU d'Esopo-Pape LLL algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    if (run_num == 0)
    {
      /* compare results */
      for (long i = 0; i < num_start_nodes; i++)
      {
        if (verify_results(distances + i*num_nodes,
                           gold_distances + i*num_nodes, num_nodes) != 0)
        {
        fprintf(stdout, "test FAILED (i = %d)\n", i);
        failed = true; anyfailed = true;
        }
      }
      if (!failed)
        fprintf(stdout, "test %d costs ok\n", num_start_nodes *num_nodes);
      
      /* check paths */
      check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                  distances, predecessors, true);
    }
  }
  fprintf(stdout,  "CPU d'Esopo-Pape LLL\talgorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_pape_lll.paths", "w");
  for (long s = 0; s < num_start_nodes; s++)
    printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */

#ifdef RUN_PAPE_LLL_SAME_MEMORY
  /* Run d'Esopo-Pape LLL algorihtm but resuse same vectors for each origin,
     like in actual use in Frank-Wolfe for TAP. Cannot do verification as
     we have done this way */
  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Running d'Esopo-Pape[reuse_same_arrays] LLL algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    sssp_pape_reuse_arrays_lll(Va, Ea, Wa, num_nodes, num_edges, start_nodes, 0,
                               distances, predecessors, num_start_nodes);

    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU d'Esopo-Pape[reuse_same_arrays] LLL algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;
  }
  fprintf(stdout,  "CPU d'Esopo-Pape[reuse_same_arrays] LLL\talgorithm average time %f ms\n", 
          total_time / NUM_RUNS);
#endif /* RUN_PAPE_LLL_SAME_MEMORY */


#ifdef RUN_SLOW_BELLMANFORD /* should implement a timeout if using this, much slower than others */
  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Running Bellman-Ford algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    for (long i = 0; i < num_start_nodes; i++)
    {
/*      fprintf(stderr, "i = %ld\n", i); */
      sssp_bellmanford(adjlist,    num_nodes, num_edges, start_nodes[i], 0,
                       distances + i*num_nodes, predecessors + i*num_nodes, 1);
    }
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU Bellman-Ford algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    if (run_num == 0)
    {
      for (long i = 0; i < num_start_nodes; i++)
      {
        if (verify_results(distances + i*num_nodes,
                           gold_distances + i*num_nodes, num_nodes) != 0)
        {
          fprintf(stdout, "test FAILED (i = %d)\n", i);
          failed = true; anyfailed = true;
        }
      }
      if (!failed)
        fprintf(stdout, "test %d costs ok\n", num_start_nodes *num_nodes);
      
      /* check paths */
      check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                  distances, predecessors, true);
    }
  }
  fprintf(stdout,  "CPU Bellman-Ford\talgorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_bellmanford.paths", "w");
  for (long s = 0; s < num_start_nodes; s++)
    printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */

#endif /*RUN_SLOW_BELLMANFORD*/


  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Running naive label-correcting algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    for (long i = 0; i < num_start_nodes; i++)
      sssp_naive_labelcorrecting(Va, Ea, Wa,  num_nodes, num_edges, start_nodes[i], 0,
                       distances + i*num_nodes, predecessors + i*num_nodes);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU naive label-correcting algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    if (run_num == 0)
    {
      for (long i = 0; i < num_start_nodes; i++)
      {
        if (verify_results(distances + i*num_nodes,
                           gold_distances + i*num_nodes, num_nodes) != 0)
        {
          fprintf(stdout, "test FAILED (i = %d)\n", i);
          failed = true; anyfailed = true;
        }
      }
      if (!failed)
        fprintf(stdout, "test %d costs ok\n", num_start_nodes *num_nodes);
      
      /* check paths */
      check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                  distances, predecessors, true);
    }
  }
  fprintf(stdout,  "CPU naive label-correcting\talgorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_naivelabelcorrecting.paths", "w");
  for (long s = 0; s < num_start_nodes; s++)
    printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */


  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Running SLF algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    for (long i = 0; i < num_start_nodes; i++)
      sssp_slf(Va, Ea, Wa, num_nodes, num_edges, start_nodes[i], 0,
                    distances + i*num_nodes, predecessors + i*num_nodes);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU SLF algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    if (run_num == 0)
    {
      for (long i = 0; i < num_start_nodes; i++)
      {
        if (verify_results(distances + i*num_nodes,
                           gold_distances + i*num_nodes, num_nodes) != 0)
        {
          fprintf(stdout, "test FAILED (i = %d)\n", i);
          failed = true; anyfailed = true;
        }
      }
      if (!failed)
        fprintf(stdout, "test %d costs ok\n", num_start_nodes *num_nodes);
      
      /* check paths */
      check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                  distances, predecessors, true);
    }
  }
  fprintf(stdout,  "CPU SLF\talgorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_slf.paths", "w");
  for (long s = 0; s < num_start_nodes; s++)
    printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */

  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Running SLF_LLL algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    for (long i = 0; i < num_start_nodes; i++)
      sssp_slf_lll(Va, Ea, Wa, num_nodes, num_edges, start_nodes[i], 0,
                    distances + i*num_nodes, predecessors + i*num_nodes);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU SLF_LLL algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    if (run_num == 0)
    {
      for (long i = 0; i < num_start_nodes; i++)
      {
        if (verify_results(distances + i*num_nodes,
                           gold_distances + i*num_nodes, num_nodes) != 0)
        {
          fprintf(stdout, "test FAILED (i = %d)\n", i);
          failed = true; anyfailed = true;
        }
      }
      if (!failed)
        fprintf(stdout, "test %d costs ok\n", num_start_nodes *num_nodes);
      
      /* check paths */
      check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                  distances, predecessors, true);
    }
  }
  fprintf(stdout,  "CPU SLF_LLL\talgorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_slf_lll.paths", "w");
  for (long s = 0; s < num_start_nodes; s++)
    printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */


#ifdef RUN_SPARSE_BLAS_BELLMAN_FORD /* FIXME there seems to be a memory leak in SparseBLAS library, don't use this */

  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Running algebraic Bellman-Ford algorithm with Sparse BLAS on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    for (long i = 0; i < num_start_nodes; i++)
    {
      bellmanford_spblas(adjlist,    num_nodes, num_edges, start_nodes[i], 
                       distances + i*num_nodes);
    }
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU sparse BLAS algebraic Bellman-Ford algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    if (run_num == 0)
    {
      for (long i = 0; i < num_start_nodes; i++)
      {
        if (verify_results(distances + i*num_nodes,
                           gold_distances + i*num_nodes, num_nodes) != 0)
        {
          fprintf(stdout, "test FAILED (i = %d)\n", i);
          failed = true; anyfailed = true;
        }
      }
      if (!failed)
        fprintf(stdout, "test %d costs ok\n", num_start_nodes *num_nodes);
      
      /* check paths */
      fprintf(stderr, "TODO no paths in algebraic Bellman-Ford yet\n");
      /* TODO implement paths in algebraic Bellman-Ford
      check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                  distances, predecessors, true);
      */
    }
  }
  fprintf(stdout,  "CPU sparse BLAS algebraic Bellman-Ford\talgorithm average time %f ms\n", 
          total_time / NUM_RUNS);


#endif /* RUN_SPARSE_BLAS_BELLMAN_FORD */

#ifdef RUN_CUSP_BELLMAN_FORD

  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Running algebraic Bellman-Ford algorithm with CUSP spmv on host CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    for (long i = 0; i < num_start_nodes; i++)
    {
      bellmanford_cusp_host(adjlist,    num_nodes, num_edges, start_nodes[i], 
                       distances + i*num_nodes);
    }
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU host CUSP sparse BLAS algebraic Bellman-Ford algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    if (run_num == 0)
    {
      for (long i = 0; i < num_start_nodes; i++)
      {
        if (verify_results(distances + i*num_nodes,
                           gold_distances + i*num_nodes, num_nodes) != 0)
        {
          fprintf(stdout, "test FAILED (i = %d)\n", i);
          failed = true; anyfailed = true;
        }
      }
      if (!failed)
        fprintf(stdout, "test %d costs ok\n", num_start_nodes *num_nodes);

      /* check paths */
      fprintf(stderr, "TODO no paths in algebraic Bellman-Ford yet\n");
      /* TODO implement paths in algebraic Bellman-Ford
      check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                  distances, predecessors, true);
      */
    }
  }
  fprintf(stdout,  "CPU host CUSP Bellman-Ford\talgorithm average time %f ms\n", 
            total_time / NUM_RUNS);

#endif /* RUN_CUSP_BELLMAN_FORD */

#endif /* ONLY_RUN_CUDA */

run_pape_cuda:
  double *pape_cuda_distances = new double[num_nodes];
  if (use_cuda)
  {
    total_time = 0.0;
    for (run_num = 0; run_num < NUM_RUNS; run_num++)
    {
      /* Run the CUDA d'Esopo Pape algorithm on all pairs */
      fprintf(stdout, "running CUDA d'Esopo-Pape algorithm on GPU...\n");
      cutilCheckError( cutResetTimer(hTimer) );
      cutilCheckError( cutStartTimer(hTimer) );

      pape_cuda(Va, Ea, Wa,  num_nodes, num_edges, start_nodes, num_start_nodes,
                 distances, predecessors);

      cutilCheckError( cutStopTimer(hTimer) );
      runtime = cutGetTimerValue(hTimer);
      fprintf(stdout,  "GPU CUDA d'Esopo-Pape et al algorithm execution time %f ms\n", runtime);
      total_time += runtime;


      /* compare results */
      failed = false;
      if (run_num == 0)
      {
        for (long i = 0; i < num_start_nodes; i++)
        {
          for (long j = 0; j < num_nodes; j++)
            pape_cuda_distances[j] = distances[j*num_start_nodes+i];
          if (verify_results(pape_cuda_distances, gold_distances + i*num_nodes,
                             num_nodes) != 0)
          {
            fprintf(stdout, "test FAILED (i = %d)\n", i);
            failed = true; anyfailed = true;
          }
        }
        if (!failed)
          fprintf(stdout, "test %d costs ok\n", num_start_nodes * num_nodes);

        /* check paths */
        check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                    distances, predecessors, false);
      }
  #ifdef TRACEBACK
      fp = fopen("test_pape_cuda.paths", "w");
      for (long s = 0; s < num_start_nodes; s++)
        printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
                   false);
      fclose(fp);
  #endif /* TRACEBACK */
    }
    fprintf(stdout,  "GPU CUDA d'Esopo-Pape\talgorithm average time %f ms\n", 
            total_time / NUM_RUNS);


    total_time = 0.0;
    for (run_num = 0; run_num < NUM_RUNS; run_num++)
    {
      fprintf(stdout, "Running algebraic Bellman-Ford algorithm with CUSP spmv on GPU...\n");
      cutilCheckError( cutCreateTimer(&hTimer) );
      cutilCheckError( cutResetTimer(hTimer) );
      cutilCheckError( cutStartTimer(hTimer) );

      for (long i = 0; i < num_start_nodes; i++)
      {
        bellmanford_cusp_device(adjlist,    num_nodes, num_edges, start_nodes[i], 
                         distances + i*num_nodes);
      }
      cutilCheckError( cutStopTimer(hTimer) );
      runtime = cutGetTimerValue(hTimer);
      fprintf(stdout,  "GPU CUSP sparse BLAS algebraic Bellman-Ford algorithm execution time %f ms\n",
              runtime);
      total_time +=  runtime;

      failed = false;
      /* compare results */
      if (run_num == 0)
      {
        for (long i = 0; i < num_start_nodes; i++)
        {
          if (verify_results(distances + i*num_nodes,
                             gold_distances + i*num_nodes, num_nodes) != 0)
          {
            fprintf(stdout, "test FAILED (i = %d)\n", i);
            failed = true; anyfailed = true;
          }
        }
        if (!failed)
          fprintf(stdout, "test %d costs ok\n", num_start_nodes *num_nodes);

        /* check paths */
        fprintf(stderr, "TODO no paths in algebraic Bellman-Ford yet\n");
        /* TODO implement paths in algebraic Bellman-Ford
        check_paths(Va,  Ea, Wa,  num_nodes, num_edges,  num_start_nodes,
                    distances, predecessors, true);
        */
      }
    }
    fprintf(stdout,  "GPU CUSP Bellman-Ford\talgorithm average time %f ms\n", 
            total_time / NUM_RUNS);

  }




  if (!anyfailed)
     printf("all tests passed\n");
  else
     printf("some test(s) failed\n");
  
  /* cleanup and exit */
  delete[] pape_cuda_distances;
  delete[] distances;
  delete[] predecessors;
  delete[] gold_distances;
  delete[] gold_predecessors;
  delete[] Va;
  delete[] Ea;
  delete[] Wa;
  cudaThreadExit();
  exit(anyfailed ? 1 : 0);
}
