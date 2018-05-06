/*****************************************************************************
 * 
 * File:    sssp_gold.cpp
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: sssp_gold.cpp 668 2011-09-08 04:40:08Z astivala $
 *
 * Get single-source shortest path on CPU for verification.
 * Uses Dijkstra's algorithm from the Boost Graph Library.
 * (Based on the example/dijkstra-example.cpp code from Boost Graph Library
 * manual
 *  http://www.boost.org/doc/libs/1_38_0/libs/graph/doc/dijkstra_shortest_paths.html)
 *
 ****************************************************************************/

#include <assert.h>

#include <cutil_inline.h>      /* CUDA SDK */

#include <boost/config.hpp>
#include <iostream>
#include <fstream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "sssp_gold.h"

using namespace boost;

/*
 * sssp_gold() -  use Dijkstra's algorithim on CPU (using Boost Graph Library)
 *                to solve single-source shortest path problem for checking.
 *
 * Parameters:
 *     adjlist - adjacnecy list (array of (node1,node2,cost) structs)
 *     num_nodes  - number of nodes 
 *     num_arcs - number of edges (entries in adjlist)
 *     source - source node 
 *     distances (OUT) - array of shortest costs from source to each node
 *     predecessors (OUT) - predecessor node on minimum spanning tree for
 *                          each node.
 * Return value:
 *     time spent running dijkstra_shortest_paths() in milliseconds
 */
double sssp_gold(adjlist_entry_t adjlist[], long num_nodes, long num_arcs,
               long source,
               double distances[], long predecessors[])
{
  unsigned int hTimer;
  typedef adjacency_list < listS, vecS, directedS,
    no_property, property < edge_weight_t, double > > graph_t;
  typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
  typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
  typedef std::pair<long, long> Edge;

  Edge *edge_array = new Edge[num_arcs];
  double *weights = new double[num_arcs];
  for (long i = 0; i < num_arcs; i++)
  {
    edge_array[i] =  Edge(adjlist[i].from, adjlist[i].to);
    weights[i] = adjlist[i].cost;
  }

  
  graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
  property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);

  std::vector<vertex_descriptor> p(num_vertices(g));
  std::vector<double> d(num_vertices(g));
  vertex_descriptor s = vertex(source, g);

  cutilCheckError( cutCreateTimer(&hTimer) );
  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );

  dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));

  cutilCheckError( cutStopTimer(hTimer) );
  double runtime = cutGetTimerValue(hTimer);

  graph_traits < graph_t >::vertex_iterator vi, vend;
  long i = 0;
  for (tie(vi, vend) = vertices(g); vi != vend; ++vi)
  {
    assert(i < num_nodes);
    distances[i] = d[*vi];
    predecessors[i] = p[*vi];
    i++;
  }

  delete[] edge_array;
  delete[] weights;

  return runtime;
}
