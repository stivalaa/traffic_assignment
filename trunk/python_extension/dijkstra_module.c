/*****************************************************************************
 * 
 * File:    dijkstra_module.c
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: dijkstra_module.c 175 2011-03-28 04:06:57Z astivala $
 *
 * Python extension to compute shortest paths using Dijkstra's algorithm
 *
 * TODO - spearetly build packed adjancany list rather than use
 * Python dictionary API in algorithm - 3 orders of magntidue overhead
 * in doing so (ie.. on ChicagoRegional network, > 2seconds to convert
 * Pyton dict do adjancency lists, but shortest paht ocmputation is only
 * a few ms)! (See pape_module.c) - using pape_module.c instead
 * of this one now.
 ****************************************************************************/

#include <Python.h>
#include <values.h> /* FLT_MAX, DBL_MAX */
#include "pqueue.h"


#define FLOATINF DBL_MAX
#define INVALID -1

/*
 * typedefs and functions for using pqueue
 */

typedef struct pqnode_t
{
        double dist; /* priority */
	long    nodenum;
	size_t pos;
        int    in_queue;
} pqnode_t;


static int
cmp_pri(double next, double curr)
{
	return (next > curr);
}


static double
get_pri(void *a)
{
	return (double) ((pqnode_t *) a)->dist;
}


static void
set_pri(void *a, double pri)
{
	((pqnode_t *) a)->dist = pri;
}


static size_t
get_pos(void *a)
{
	return ((pqnode_t *) a)->pos;
}


static void
set_pos(void *a, size_t pos)
{
	((pqnode_t *) a)->pos = pos;
}
	

/*
 * shortest path function
 */

static PyObject* dijkstra(PyObject* self, PyObject* args)
{
  PyObject *graph_dict;
  long start_node;
  PyObject *node_dict;
  long *prev;
  long num_nodes;
  long v;
  pqnode_t *unode;
  pqueue_t *queue;
  pqnode_t *qnode;
  PyObject *neighbour_list;
  Py_ssize_t i;
  PyObject *vkey, *uvdist_obj;
  double uvdist, newcost;
  PyObject *pred_list;

  if (!PyArg_ParseTuple(args, "Oll", &graph_dict, &num_nodes,
                        &start_node))
    return NULL;

  if (start_node >= num_nodes)
  {
    PyErr_SetString(PyExc_IndexError, "start_node >= num_nodes");
    return NULL;
  }

  if (!(prev = (long *)malloc(num_nodes * sizeof(long))))
    return PyErr_NoMemory();
  
  if (!(qnode = (pqnode_t *)malloc(num_nodes *sizeof(pqnode_t))))
     return PyErr_NoMemory();
  
  queue = pqueue_init(num_nodes, cmp_pri, get_pri, set_pri, get_pos, set_pos);
  if (!queue)
    return PyErr_NoMemory();
    

  for (v = 0; v < num_nodes; v++)
  {
    prev[v] = INVALID;
    qnode[v].nodenum = v;
    qnode[v].in_queue = 1;
    if (v == start_node)
      qnode[v].dist = 0;
    else
      qnode[v].dist = FLOATINF;
    pqueue_insert(queue, &qnode[v]);
  }

  while (pqueue_size(queue) > 0)
  {
    unode = (pqnode_t *)pqueue_pop(queue);  /* node with smallest dist */
    unode->in_queue = 0;
    if (unode->dist == FLOATINF)
      break;
    /* get neighbours of unode */
    node_dict = PyDict_GetItem(graph_dict, PyInt_FromLong(unode->nodenum));
    if (!node_dict)
      continue; /* no neighbours */
    neighbour_list = PyDict_Keys(node_dict);
    for (i = 0; i < PyList_GET_SIZE(neighbour_list); i++)
    {
      vkey = PyList_GET_ITEM(neighbour_list, i);
      v = PyInt_AsLong(vkey);
      if (!qnode[v].in_queue)
        continue;
      uvdist_obj = PyDict_GetItem(node_dict, vkey);
      if (!uvdist_obj) /* should never happen */
      {
        PyErr_Format(PyExc_KeyError, "graph[%ld][%ld] no value\n",
                        unode->nodenum, v);
        return NULL;
      }
      uvdist = PyFloat_AsDouble(uvdist_obj);
      newcost = unode->dist + uvdist;
/*      printf("FIXME u = %ld v = %ld uvdist = %f newcost = %f\n", unode->nodenum, v, uvdist, newcost); */
      if (newcost < qnode[v].dist)
      {
        pqueue_change_priority(queue, newcost, &qnode[v]);
        prev[v] = unode->nodenum;
      }
    }
  }

  pred_list = PyList_New(0);
  for (v = 0; v < num_nodes; v++)
  {
    PyList_Append(pred_list, PyInt_FromLong(prev[v]));
  }

  pqueue_free(queue);
  free(prev);
  free(qnode);

  return pred_list;
}

static PyMethodDef DijkstraMethods[] =
{
  {"dijkstra", dijkstra, METH_VARARGS, 
   "Single-source shortest path.\n\n"
   "shortest_path(graph, num_nodes, start_node)\n"
   "\n"
   "Single-source shortest path from start_node in graph,\n"
   "using Dijkstra' algorithm.\n"
   "Parameters:\n"
   "  graph - dict of dicts where G[v][w] for any v,w is cost v->w\n"
   "          v and w must be integers, cost must be double\n"
   "  num_nodes - number of nodes in graph\n"
   "  start_node int for start to find path\n"
   "Return value:\n"
   "  list pred[i] where pred[i] is predecessor of node i on shortest path\n"
  },
  {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
initcdijkstra(void)
{
  (void) Py_InitModule("cdijkstra", DijkstraMethods);
}
