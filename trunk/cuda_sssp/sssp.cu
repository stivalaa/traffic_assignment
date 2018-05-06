/*****************************************************************************
 * 
 * File:    sssp.cu
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: sssp.cu 260 2011-04-28 04:29:24Z astivala $
 *
 * Single-source shortest path implementation using CUDA based on
 * 
 * Okuyama, Ino, Hagihara 2008 "A Task Parallel Algorithm for Computing
 * the costs of All-Pairs Shortest Paths on the CUDA-compatible GPU"
 * Intl. Symp. Parallel Distributed Processing with Applications. 284-291
 *
 * See also:
 *
 * Harish and Narayanan 2007 "Accelerating Large Graph Algorithms on the GPU
 * Using CUDA" HiPC 2007, LNCS 4873: 197-208
 * 
 * and:
 * 
 * Martin, Torres, Gavilanes "CUDA Solutions for the SSSP Problem"
 * ICCS 2009, LNCS 5544: 904-913
 *
 ****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <cutil_inline.h>      /* CUDA SDK */

#include "sssp.h"
#include "sssp_gold.h"
#include "sssp_pape.h"
#include "harish_host.h"
#include "okuyama_host.h"
#include "pape_cuda_host.h"

#undef DUMP_RESULTS
#undef TRACEBACK

const float EPS = 1e-08;  // error toelrance

/*****************************************************************************
 *
 * Local functions
 *
 ****************************************************************************/

//Align a to nearest higher multiple of b
extern "C" int iAlignUp(int a, int b){
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
void dump_packed_arrays(const int Va[], const int Ea[], const float Wa[],
                        int num_nodes, int num_edges)
{
  int i;
  for (i = 0; i < num_nodes; i++)
    printf("%d ", Va[i]);
  printf("\n");
  for (i = 0; i < num_edges; i++)      
    printf("%4d ", Ea[i]);
  printf("\n");
  for (i = 0; i < num_edges; i++)
    printf("%4.2f ", Wa[i]);
  printf("\n");
}
#endif  /* DEBUG */


/*
 * verify_results() -compare two floating point arrays 
 *
 * Parameters:
 *    test -  results to check
 *    ref - reference versino to compare against 
 *    n - length of arrays to compare
 *
 * Return value:
 *    0 if OK else -1
 */
static int verify_results(float *test, float *ref, int n)
{
  int rc = 0;
  for (int i = 0; i < n; i++)
    if (fabsf(ref[i] - test[i]) > EPS)
    {
      fprintf(stdout, "i = %d ref = %f test = %f\n", i, ref[i], test[i]);
      rc = -1;
    }
  return rc;
}



/* Retreive the mincost path from i to j implied by the
   distance and predecessor matrices
   NB this path is stored backwards (from j to i) */
int Min_getpath(float *distances, int *predecessors, int num_nodes,
                int num_start_nodes, int i, int j, int **path, 
                int *pathlen, bool sourcemajor)
{
  int idx = sourcemajor ? i * num_nodes + j : j * num_start_nodes + i;
  if (distances[idx] == FLOATINF)
    return 0; // no path
  int k = j;
  while (k != i)
  {
    assert(*pathlen < num_nodes);
    **path = k;
    (*path)++;
    (*pathlen)++;
    int k_idx = sourcemajor ? i * num_nodes + k : k * num_start_nodes + i;
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
static int verify_path(int path[], int pathlen, float cost,
                       int Va[], int Ea[], float Wa[], int num_nodes,
                       int num_edges)
{
  if (cost == FLOATINF && pathlen < 1)
    return 0; //no path
  float acost = 0;
  for (int k = pathlen-1; k > 0; k--)
  {
    if (path[k] > num_nodes)
    {
      fprintf(stderr, "bad path: %d > num_nodes at %d\n", path[k], k);
      return -1;
    }
    if (path[k-1] > num_nodes)
    {
      fprintf(stderr, "bad path: %d > num_nodes at %d\n", path[k+1], k+1);
      return -1;
    }

    int start = Va[path[k]];
    int end = Va[path[k]+1];
    int i;
    for (i = start; i < end && Ea[i] != path[k-1]; i++)
      /*nothing*/;
    if (Ea[i] != path[k-1])
    {
      fprintf(stderr, "bad path: no edge from %d to %d\n", path[k],
              path[k-1]);
      return -1;
    }
    acost += Wa[i];
  }
  if (fabsf(acost - cost) > EPS)
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
static int check_paths(int Va[], int Ea[], float Wa[], int num_nodes,
                       int num_edges, int num_start_nodes,
                       float distances[],
                       int predecessors[], bool sourcemajor)
{

  /* check paths */
  bool pathtest_failed = false;
  int *thepath  = new int[num_nodes];
  for (int i = 0; i < num_start_nodes; i++)
  {
    for (int j = 0; j < num_nodes; j++)
    {
      int *path = thepath;
      int **pathptr = &path;
      int pathlen = 0;
      Min_getpath(distances, predecessors, num_nodes, num_start_nodes,
                  i, j, pathptr, &pathlen, sourcemajor);
      int idx = sourcemajor ? i * num_nodes + j : j * num_start_nodes + i;
      if (verify_path(thepath, pathlen, 
                      distances[idx],
                      Va, Ea, Wa, num_nodes,
                      num_edges) != 0)
      {
        printf("verify path FAILED source %d dest %d\n", i, j);
        pathtest_failed = true;
      }
    }
  }
  if (!pathtest_failed)
    printf("verify %d paths ok\n", num_start_nodes*num_nodes);
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
static int read_adjlist(FILE *fp, adjlist_entry_t **adjlist, int *num_nodes)
{
  int x,y,nnz,n;
  adjlist_entry_t *entry;
  int max_node_num = 0;

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
    if (fscanf(fp, "%d\t%d\t%f\n",&entry->from, &entry->to, &entry->cost)!= 3)
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
      max_node_num = max(entry->from, entry->to);

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
                                     int num_nodes,
                                     int num_edges,
                                     int Va[],
                                     int Ea[],
                                     float Wa[])
{
  // sort by 'from' node ascending and within that by 'to' node ascending
  qsort(adjlist, num_edges, sizeof(adjlist_entry_t), adjlist_entry_compar);
  int v = 0;  // index into Va
  int e = 0;  // index into Ea and Wa and adjlist

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
 * Return value: None
 *
 * Sets the global use_shared_memory and calls cudaSetDevice()
 * to set the GPU to use.
 */
static void choose_gpu(int major_req, int minor_req)
{
/*
    int devnum = cutGetMaxGflopsDeviceId();
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
//        use_shared_memory = false;
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
//          use_shared_memory = true;
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
static void printpaths(FILE *fp, float cost[], int predecessor[], int source,
                       int num_nodes, int num_start_nodes, bool sourcemajor)
{
  for (int dest = 0; dest < num_nodes; dest++)
  {
    int idx = sourcemajor ? source*num_nodes+dest : dest*num_start_nodes+source;
      
    if (cost[idx] == FLOATINF)
    {
      fprintf(fp, "%d -> %d: no path\n", source+1, dest+1);
      continue;
    }
    fprintf(fp, "%d -> %d (%f): ", source+1, dest+1, cost[idx]);
    int k = dest;
    int k_idx = sourcemajor ? source*num_nodes+k : k*num_start_nodes+source;
    while (k != source)
    {
      fprintf(fp, "%d <- ", k+1);
      k = predecessor[k_idx];
      k_idx = sourcemajor ? source*num_nodes+k : k*num_start_nodes+source;
    }
    fprintf(fp, "%d \n", source+1);
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
  float *distances;
  int *predecessors;
  bool failed = false, anyfailed = false;
  const int NUM_RUNS = 10; // run each test this many times to avg timings


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
  adjlist_entry_t *adjlist; int num_nodes = 0;
  int  num_edges = read_adjlist(fp, &adjlist, &num_nodes);
  if (num_edges < 0)
  {
    fprintf(stderr, "reading input file %s failed\n", inputfilename);
    exit(1);
  }

  /* Initialize GPU */
  // need at least compute capability 1.2 for 64 bit atomic CAS instruction
  choose_gpu(1,2);  /* calls cudaSetDevice() */

  /* Run the CUDA SSSSP imlementation */

  /* first we have to convert adjacnecy lists to packed array format */

  // pad the arrays out to a multiple of BLOCK_SIZE
//  int alloc_num_nodes =  iAlignUp(num_nodes, BLOCK_SIZE);
  int alloc_num_nodes = num_nodes+1; //no longer need the above
  int *Va = new int[alloc_num_nodes];
  int *Ea = new int[num_edges];
  float *Wa = new float[num_edges];
//  for (int i = num_nodes; i < alloc_num_nodes; i++)
//    Va[i] = Va[num_nodes]+1;  // no edges on pad nodes

  adjlist_to_packed_arrays(adjlist, num_nodes, num_edges, Va, Ea, Wa);

#ifdef DEBUG
  dump_packed_arrays(Va, Ea, Wa, num_nodes, num_edges);
#endif  

//  const int NUM_START_NODES =  147; // # of zones in winnipeg data // num_nodes; // / 2; //; // XXX
  const int NUM_START_NODES =  1790; // XXX # of zones in ChicagoRegional
//  const int NUM_START_NODES = num_nodes / 2;

  int num_start_nodes = NUM_START_NODES;

  assert(num_start_nodes <= num_nodes);

  int *start_nodes = new int[num_start_nodes];
  for (int i = 0; i < num_start_nodes; i++)
    start_nodes[i]  = i;

  float *gold_distances = new float[num_nodes*num_start_nodes];
  int *gold_predecessors = new int[num_nodes*num_start_nodes];


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
    for (int i = 0; i < num_start_nodes; i++)
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
  fprintf(stdout,  "CPU Dijkstra algorithm average time %f ms\n", 
          total_time / NUM_RUNS);
  fprintf(stdout,  "CPU Dijkstra algorithm average internal time %f ms\n", 
          total_internal_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("reference.paths", "w");
  for (int s = 0; s < num_start_nodes; s++)
    printpaths(fp, gold_distances, gold_predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */


  distances = new float[num_nodes*num_start_nodes];
  predecessors = new int[num_nodes*num_start_nodes];

  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Runing d'Esopo-Pape algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    for (int i = 0; i < num_start_nodes; i++)
      sssp_pape(Va, Ea, Wa, num_nodes, num_edges, start_nodes[i], 0,
                distances + i*num_nodes, predecessors + i*num_nodes);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU d'Esopo-Pape algorithm execution time %f ms\n", runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    for (int i = 0; i < num_start_nodes; i++)
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
  fprintf(stdout,  "CPU d'Esopo-Pape algorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_pape.paths", "w");
  for (int s = 0; s < num_start_nodes; s++)
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

    for (int i = 0; i < num_start_nodes; i++)
      sssp_pape_lll(Va, Ea, Wa, num_nodes, num_edges, start_nodes[i], 0,
                    distances + i*num_nodes, predecessors + i*num_nodes);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU d'Esopo-Pape LLL algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    for (int i = 0; i < num_start_nodes; i++)
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
  fprintf(stdout,  "CPU d'Esopo-Pape LLL algorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_pape_lll.paths", "w");
  for (int s = 0; s < num_start_nodes; s++)
    printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */



  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "Running Bellman-Ford algorithm on CPU...\n");
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    for (int i = 0; i < num_start_nodes; i++)
      sssp_bellmanford(Va, Ea, Wa, num_nodes, num_edges, start_nodes[i], 0,
                    distances + i*num_nodes, predecessors + i*num_nodes);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU Bellman-Ford algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    for (int i = 0; i < num_start_nodes; i++)
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
  fprintf(stdout,  "CPU Bellman-Ford algorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_bellmanford.paths", "w");
  for (int s = 0; s < num_start_nodes; s++)
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

    for (int i = 0; i < num_start_nodes; i++)
      sssp_slf(Va, Ea, Wa, num_nodes, num_edges, start_nodes[i], 0,
                    distances + i*num_nodes, predecessors + i*num_nodes);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU SLF algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    for (int i = 0; i < num_start_nodes; i++)
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
  fprintf(stdout,  "CPU SLF algorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_slf.paths", "w");
  for (int s = 0; s < num_start_nodes; s++)
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

    for (int i = 0; i < num_start_nodes; i++)
      sssp_slf_lll(Va, Ea, Wa, num_nodes, num_edges, start_nodes[i], 0,
                    distances + i*num_nodes, predecessors + i*num_nodes);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "CPU SLF_LLL algorithm execution time %f ms\n",
            runtime);
    total_time +=  runtime;

    failed = false;
    /* compare results */
    for (int i = 0; i < num_start_nodes; i++)
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
  fprintf(stdout,  "CPU SLF_LLL algorithm average time %f ms\n", 
          total_time / NUM_RUNS);

#ifdef TRACEBACK
  fp = fopen("test_slf_lll.paths", "w");
  for (int s = 0; s < num_start_nodes; s++)
    printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
               true);
  fclose(fp);
#endif /* TRACEBACK */



  float *okuyama_distances = new float[num_nodes];

  total_time = 0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    fprintf(stdout, "running Harish & Narayanan algorithm on GPU...\n");
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );
    
    /* Run the Harish & Narayanan algorithm */
    for (int i = 0; i < num_start_nodes; i++)
      harish_sssp(Va, Ea, Wa, alloc_num_nodes, num_nodes, num_edges,
                  start_nodes[i], distances + i*num_nodes,
                  predecessors + i*num_nodes); 
    
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "GPU Harish&Narayanan algorithm execution time %f ms\n", runtime);
    total_time += runtime;
    

    failed = false;
    /* compare results */
    for (int i = 0; i < num_start_nodes; i++)
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

#ifdef DUMP_RESULTS
    /* dump distances and predcessors */
    printf("ref:\n");
    for (int i = 0; i < num_nodes; i++)
    {
    if (gold_predecessors[i] != i) // if ==, either source or unreachable
    {
      printf("%d: %d (%f)\n", i+1, gold_predecessors[i]+1, gold_distances[i]);
    }
    }
    printf("test:\n");
    
#ifdef DEBUG
    for (int i  =0; i < num_start_nodes; i++)
    {
      for (int j =0; j < num_nodes; j++)
      {
        printf("%4d ", predecessors[i*num_start_nodes+j]);
      }
      printf("\n");
    }
#endif /* DEBUG */      
    for (int i = 0; i < num_nodes; i++)
    {
      if (predecessors[i] != i) // if ==, either source or unreachable
      {
        printf("%d: %d (%f)\n", i+1, predecessors[i]+1, distances[i]);
      }
    }
#endif /* DUMP_RESULTS */
    
  } 
  fprintf(stdout,  "GPU Harish&Narayanan algorithm average time %f ms\n",
          total_time / NUM_RUNS);
 
  memset(distances, 0, num_nodes*num_start_nodes*sizeof(distances[0]));
  memset(predecessors, 0, num_nodes*num_start_nodes*sizeof(predecessors[0]));

run_okuyama:
  total_time = 0.0;
  for (run_num = 0; run_num < NUM_RUNS; run_num++)
  {
    /* Run the Okuyama et al algorithm on all pairs */
    fprintf(stdout, "running Okuyama et al algorithm on GPU...\n");
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );
    
    okuyama_sssp(Va, Ea, Wa,  num_nodes, num_edges, start_nodes, num_start_nodes,
               distances, predecessors);
    
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stdout,  "GPU Okuyama et al algorithm execution time %f ms\n", runtime);
    total_time += runtime;
    
    /* compare results */
    failed = false;
    for (int i = 0; i < num_start_nodes; i++)
    {
      for (int j = 0; j < num_nodes; j++)
        okuyama_distances[j] = distances[j*num_start_nodes+i];
      if (verify_results(okuyama_distances, gold_distances + i*num_nodes,
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
    
#ifdef TRACEBACK
    fp = fopen("test.paths", "w");
    for (int s = 0; s < num_start_nodes; s++)
      printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
                 false);
    fclose(fp);
#endif /* TRACEBACK */
  }
    fprintf(stdout,  "GPU Okuyama et al algorithm average time %f ms\n", 
            total_time / NUM_RUNS);


run_pape_cuda:
  float *pape_cuda_distances = new float[num_nodes];
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
    for (int i = 0; i < num_start_nodes; i++)
    {
      for (int j = 0; j < num_nodes; j++)
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
    
#ifdef TRACEBACK
    fp = fopen("test_pape_cuda.paths", "w");
    for (int s = 0; s < num_start_nodes; s++)
      printpaths(fp, distances, predecessors, s, num_nodes, num_start_nodes,
                 false);
    fclose(fp);
#endif /* TRACEBACK */
  }
  fprintf(stdout,  "GPU CUDA d'Esopo-Pape algorithm average time %f ms\n", 
          total_time / NUM_RUNS);


  if (!anyfailed)
     printf("all tests passed\n");
  else
     printf("some test(s) failed\n");
  
  /* cleanup and exit */
  delete[] pape_cuda_distances;
  delete[] okuyama_distances;
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
