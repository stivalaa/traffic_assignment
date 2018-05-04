/*****************************************************************************
 * 
 * File:    pape_module.c
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: pape_module.c 784 2011-10-05 04:31:45Z astivala $
 *
 * Python extension to compute single-source shortest paths using 
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
 ****************************************************************************/

#include <Python.h>
#include <values.h> /* FLT_MAX, DBL_MAX */
#include <sys/time.h>
#include <sys/resource.h>

#undef DEBUG
#undef TIME_DEBUG

/****************************************************************************
 *
 * constants and type definitions
 *
 ****************************************************************************/

#define FLOATINF      DBL_MAX
#define INVALID      -1
#define WAS_IN_QUEUE -2

/* adjacency list entry */
typedef struct adjlist_entry_s 
{
    int from;          /* from node id (0..n-1) */
    int to;            /* to node id (0..n-1) */
    double cost;        /* cost of arc from 'from' to 'to' >= 0 */
} adjlist_entry_t;



/****************************************************************************
 *
 * constants and type definitions
 *
 ****************************************************************************/

/* DODGY: we will make these static so pape.build_adjlist() can be called
   to build them once and weights update with pape.update_weights()
   and the actual shortgst path done with pape.pape() not actually
   passing graph to the latter at all (using these static arrays instead(
   TODO:  make an opaque handle or something (too hard to bother) 

*/
static long *Va, *Ea;
static double *Wa;
static long num_nodes;
static long global_num_edges;
static long first_thru_node;

/****************************************************************************
 *
 * local functions to convert python dictionary to packed adjacency list
 *
 ****************************************************************************/
	

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
 * The adjlist input is sorted in place.
 */
static void adjlist_to_packed_arrays(adjlist_entry_t adjlist[],
                                     long num_nodes,
                                     long num_edges,
                                     long Va[],
                                     long Ea[],
                                     double Wa[])
{
  /* sort by 'from' node ascending and within that by 'to' node ascending */
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

/*
 * build_packed_adjlist_from_dict() - build pakced adjacney list from Python dict
 *
 * Parameters:
 *        graph_dict - graph as Python dictionary graph_dict[u][v] = cost
 *        num_nodes - number of nodes in graph (>= largest node# +1)
 *        Va[] (OUT)-array of indices to head of each adj list in Ea 
 *        Ea[] (OUT)-each entry is 'to' node in list for 'from' node
 *        Wa[] (OUT)-each entry is cost of corresponding Ea entry
 * 
 * Return value:
 *       0 if Ok else nonzero
 *
 *
 * The adjlist input must be sorted by from and then to node ascending
 * on input.
 */
static int build_packed_adjlist_from_dict(PyObject *graph_dict,
                                          long num_nodes,
                                          long **Va, long **Ea, double **Wa)
{
  const int INITIAL_NUM = 1000; /* allocate this many to begin with */
  long num_edges;
  PyObject *n_list;
  adjlist_entry_t *adjlist;
  PyObject *neighbour_list;
  PyObject *vkey, *uvdist_obj, *ukey;
  PyObject *node_dict;
  double uvdist;
  Py_ssize_t i,j;
  long u,v;
#ifdef TIME_DEBUG
  struct rusage starttime,endtime;
  int otime;
#endif

#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  if (!(adjlist = (adjlist_entry_t *)malloc(INITIAL_NUM*sizeof(adjlist_entry_t))))
    return -1;
  num_edges = 0;
  n_list = PyDict_Keys(graph_dict); /* list of all nodes with out edges */
  for (i = 0; i < PyList_GET_SIZE(n_list); i++)
  {
    ukey = PyList_GET_ITEM(n_list, i);
    u = PyInt_AsLong(ukey);
    /* get neighbours of u */
    node_dict = PyDict_GetItem(graph_dict, ukey);
    if (!node_dict)
      continue; /* no neighbours */
    neighbour_list = PyDict_Keys(node_dict);
    for (j = 0; j < PyList_GET_SIZE(neighbour_list); j++)
    {
      vkey = PyList_GET_ITEM(neighbour_list, j);
      v = PyInt_AsLong(vkey);
      uvdist_obj = PyDict_GetItem(node_dict, vkey);
      if (!uvdist_obj) /* should never happen */
      {
        PyErr_Format(PyExc_KeyError, "graph[%ld][%ld] no value\n",  u, v);
        return -1;
      }
      uvdist = PyFloat_AsDouble(uvdist_obj);
/*      printf("u = %ld v = %ld uvdist = %f\n", u, v, uvdist); */
      if (num_edges >= INITIAL_NUM)
        if (!(adjlist = (adjlist_entry_t *)realloc(adjlist, 
                                     (num_edges + 1)*sizeof(adjlist_entry_t))))
          return -1;
      adjlist[num_edges].from = u;
      adjlist[num_edges].to = v;
      adjlist[num_edges].cost = uvdist;
      num_edges++;
    }
  }
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &endtime);
  otime = 1000 * endtime.ru_utime.tv_sec + endtime.ru_utime.tv_usec/1000
          + 1000 * endtime.ru_stime.tv_sec + endtime.ru_stime.tv_usec/1000
          - 1000 * starttime.ru_utime.tv_sec - starttime.ru_utime.tv_usec/1000
          - 1000 * starttime.ru_stime.tv_sec - starttime.ru_stime.tv_usec/1000;
  printf("adjlist build time: %d ms\n", otime);
#endif

  if (!(*Va = (long *)malloc((num_nodes+1) * sizeof(long))))
    return -1;
  if (!(*Ea = (long *)malloc(num_edges * sizeof(long))))
    return -1;
  if (!(*Wa = (double *)malloc(num_edges * sizeof(double))))
    return -1;

#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &starttime);
#endif
  adjlist_to_packed_arrays(adjlist, num_nodes, num_edges, *Va, *Ea, *Wa);
#ifdef TIME_DEBUG
  getrusage(RUSAGE_SELF, &endtime);
  otime = 1000 * endtime.ru_utime.tv_sec + endtime.ru_utime.tv_usec/1000
          + 1000 * endtime.ru_stime.tv_sec + endtime.ru_stime.tv_usec/1000
          - 1000 * starttime.ru_utime.tv_sec - starttime.ru_utime.tv_usec/1000
          - 1000 * starttime.ru_stime.tv_sec - starttime.ru_stime.tv_usec/1000;
  printf("packed adjlist build time: %d ms\n", otime);
#endif
  free(adjlist); adjlist = NULL; /* no longer need this, only packed arrays */
  
#ifdef DEBUG
  dump_packed_arrays(*Va, *Ea, *Wa, num_nodes, num_edges);
#endif  

  global_num_edges = num_edges;

  return 0;
}

/****************************************************************************
 *
 * single source shortest path function callable from python
 *
 ****************************************************************************/

/*
 * build internal packed adjlist from python dict
 */
static PyObject* build_adjlist(PyObject *self, PyObject *args)
{
  PyObject *graph_dict;

  if (!PyArg_ParseTuple(args, "Oll", &graph_dict,  &num_nodes, &first_thru_node))
    return NULL;

  /* convert python dictionary representing graph to packed adjacency list */
  if (build_packed_adjlist_from_dict(graph_dict, num_nodes, &Va, &Ea, &Wa) != 0)
    return PyErr_NoMemory();

/*  printf("build_adjlist num_nodes = %ld\n", num_nodes); /* XXX */

  Py_RETURN_NONE;
}


/*
 * update weights in internal packed adjlist
 */
static PyObject* update_weights(PyObject *self, PyObject *args)
{
  PyObject *weight_list;
  Py_ssize_t i;

  if (!PyArg_ParseTuple(args, "O", &weight_list))
    return NULL;
  
  if (PyList_GET_SIZE(weight_list) != global_num_edges)
  {
    PyErr_SetString(PyExc_IndexError, "list length != num edges");
    return NULL;
  }

  for (i = 0; i < PyList_GET_SIZE(weight_list); i++)
    Wa[i] =  PyFloat_AsDouble(PyList_GET_ITEM(weight_list, i));

  Py_RETURN_NONE;
}

/*
 * run d'Esopo-Pape algorithm to get shortest paths 
 */
static PyObject* pape(PyObject* self, PyObject* args)
{

  long start_node;
  long *prev;
  long u,v;
  long i;
  long *queue_next;
  long queue_first, queue_last;
  double uvdist, newcost;
  PyObject *pred_list;
  double *dist;
  struct rusage starttime,endtime;
  int otime;

  if (!PyArg_ParseTuple(args, "l",  &start_node))
    return NULL;

/*  printf("pape num_nodes = %ld\n", num_nodes); /* XXX */


  if (start_node >= num_nodes)
  {
    PyErr_SetString(PyExc_IndexError, "start_node >= num_nodes");
    return NULL;
  }

  if (!(prev = (long *)malloc(num_nodes * sizeof(long))))
    return PyErr_NoMemory();

  if (!(dist = (double *)malloc(num_nodes * sizeof(double))))
    return PyErr_NoMemory();

  if (!(queue_next = (long *)malloc(num_nodes *sizeof(long))))
     return PyErr_NoMemory();
  
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
      uvdist = Wa[i]; /* with this cost on edge */
      newcost = dist[u] + uvdist;
      if (newcost < dist[v])
      {
        dist[v] = newcost;
        prev[v] = u;
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
  pred_list = PyList_New(0);
  for (v = 0; v < num_nodes; v++)
  {
    PyList_Append(pred_list, PyInt_FromLong(prev[v]));
  }

  free(queue_next);
  free(prev);
  free(dist);

   /* FIXME there is a memory leak in this fuction (NOT just the one
    *  below with commented out free(), some probably with refcounting
    *  meaning the pred_list of intermediate PyInt objects or something
    *  are building up and never being deleted... 
    */

  /* TODO free memory somewhere when finished */
  /* free(Va); */
  /* free(Ea); */
  /* free(Wa); */ 

  return pred_list;
}

static PyMethodDef PapeMethods[] =
{
  {"pape", pape, METH_VARARGS, 
   "Single-source shortest path.\n\n"
   "pape(graph, num_nodes, start_node)\n"
   "\n"
   "Single-source shortest path from start_node in graph,\n"
   "using d'Esopo-Pape algorithm.\n"
   "Parameters:\n"
/*   "  graph - dict of dicts where G[v][w] for any v,w is cost v->w\n"
   "          v and w must be integers, cost must be double\n"
*/
/*   "  num_nodes - number of nodes in graph\n" */
   "  start_node int for start to find path\n"
   "Return value:\n"
   "  list pred[i] where pred[i] is predecessor of node i on shortest path\n"
  },
  {"build_adjlist", build_adjlist, METH_VARARGS,
   "Build internal packed adjacency list from dict\n"
   "Prameters:\n"
   "  graph - dict of dicts where G[v][w] for any v,w is cost v->w\n"
   "          v and w must be integers, cost must be double\n"
   "  num_nodes - number of nodes in graph\n"
   "  first_thru_node - first node that paths allowed to pass through (apart from start nodes)\n"
  },
  {"update_weights", update_weights, METH_VARARGS,
   "Update weights of internal  packed adjacney list\n"
   "Parameters:\n"
   "  weights - list of weight values for each edge\n"
   "            sorted so that edges are oredered by from_node ascending\n"
   "            and within that to_node ascending\n"
  },
  {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
initpape(void)
{
  (void) Py_InitModule("pape", PapeMethods);
}
