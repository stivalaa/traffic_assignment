/*****************************************************************************
 * 
 * File:    tap_functions.c
 * Author:  Alex Stivala
 * Created: February 2011
 *
 * $Id: tap_functions.c 525 2011-08-08 04:35:34Z astivala $
 *
 * Traffic assignment functions
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
#include "parsetapfiles.h"
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


/*
 * bpr_cost_function()
 *
 * BPR function: travel_time(Q) = T_0 * (1 + alpha*(Q/Q_max)^beta) 
 *
 * Parameters:
 *      link  - link parameter struct
 *      Q     - volume on link
 *  
 * Return value:
 *      cost on link given volume (Q) using the link's parameters
 *      in the Bureau of Public Roads (BPR) function
 */
double bpr_cost_function(const link_data_t *link, double Q)
{
  double alpha = link->B;
  double beta  = link->power;
  double Q_max = link->capacity;
  assert(Q_max > 0);
  return link->free_flow_time * (1.0 + alpha * pow(Q / Q_max, beta));
  /* return link->free_flow_time * (1.0 + alpha *  */
  /*                                (link->B / pow(link->capacity, link->power)) *  */
  /*                                pow(Q, beta)); */
}



/*
 * link_directionral_derivative() - directional derivative of objective function
 *
 * Used for line search in Frank-Wolfe algorithm.
 *
 * Parameters:
 *      net - net structure parsed from net file
 *      volume[] - link volumes
 *      delta_volume - auxilliary volumes (delta)
 *      lambda - step size
 *
 * Return value:
 *    directional derivative of objective function for link-based UE
 */
double link_directional_derivative(const net_data_t *net,
                                   double volume[],
                                   double delta_volume[], double lambda)
{
  long k;
  double link_cost_sum = 0;
  double vol;

  for (k = 0; k < net->num_links; k++)
  {
    vol = volume[k] + lambda * delta_volume[k];
    link_cost_sum += bpr_cost_function(&net->links[k], vol) * delta_volume[k];
  }
  return link_cost_sum;
}

/*
 * linkcost_integral() - integral of link cost function on a single link
 *
 * Used to compute the objective function for UE with Frank-Wolfe
 * with BPR cost function
 *
 * Parameters:
 *      link  - link parameter struct
 *      Q     - volume on link
 *    
 * Return value:
 *     value of integral of link cost on this link
 */
double linkcost_integral(const link_data_t *link, double Q)
{
  double integral;
  if (link->power >= 0)
  {
    integral =  Q * link->free_flow_time *
      ( 1 + ( (link->B / pow(link->capacity, link->power))  / (link->power+1) )
        * pow(Q, link->power) );
  }
  else
  {
    integral = 0;
  }
  return integral;
}


/*
 * links_objective_function() - objective function for UE
 *
 * Used to compute the objective function for UE with Frank-Wolfe
 * with BPR cost function
 *
 *
 * Parameters:
 *      net - net structure parsed from net file
 *      volume[] - link volumes
 *
 * Return value:
 *    sum of link cost integrals 
 */
double links_objective_function(const net_data_t *net, double volumes[])
{
  long k;
  double total = 0;

  for (k = 0; k < net->num_links; k++)
    total += linkcost_integral(&net->links[k], volumes[k]) ;
  
  return total;
}

/*
 * update_link_costs()
 *
 * update the cost of each link using the BPR cost function and volume
 *
 * Parameters:
 *    net - net structure parsed from net file
 *    link_volumes - volume on each link
 *    link_costs  (OUT) - cost on each link
 *
 * Return value:
 *    None.
 */
void update_link_costs(net_data_t *net, double link_volumes[],
                       double link_costs[])
{
  int k;
  for (k = 0; k < net->num_links; k++)
  {
    link_costs[k] = bpr_cost_function(&net->links[k], link_volumes[k]);
  }
}

/*
 * cost_step()
 *
 * Given link struct and volume for that link, return the step (in the 
 * division of volume up to capacity that the volume is at step
 * function at that volume
 *
 * Parameters:
 *    link - link structure
 *    vol  - volume on the link
 *
 * Return value:
 *    step (0, 1, .. NUM_STEPS) that the volume is at on that link
 */
double cost_step(link_data_t *link, double vol)
{
  double flowstepsize;
  double step;
  if (vol > link->capacity)
    return NUM_STEPS;
  flowstepsize = link->capacity / (double)NUM_STEPS;
  step = floor(vol / flowstepsize); 
  return step;
}


/*
 * cost_step_function()
 *
 * Given a link struct and volume on that link, return the value of the
 * cost step function at that volume
 * 
 * Parameters;
 *    link - link struct
 *    vol - current volume on link
 *
 * Return value
 *    cost at that volume on link
 */
double cost_step_function(link_data_t *link, double vol)
{
  double step;
  double cost;
  double flowstepsize = link->capacity / (double)NUM_STEPS;
  step = cost_step(link, vol);
  if (vol > link->capacity)
    cost = bpr_cost_function(link, link->capacity);
  else
    cost = bpr_cost_function(link, step * flowstepsize);
  return cost;
}


/* 
 * distance_from_next_step()
 *
 * Return the distance (volume difference) from current volume on link
 * to the point on a link where it will go up to the next step in the
 * cost step function.
 *
 * Parameters:
 *    links - link data array
 *    link_volumes - link volume array
 *    link_idx - index of link (net.links[link_idx], link_volumes[link_idx])
 *    
 * Return value:
 *    distance (additional volume required) to get to nest step in cost step fn
 */
double distance_from_next_step(link_data_t links[], 
                                      double link_volumes[],
                                      long link_idx)
{
  double flowstepsize;
  double next_step_vol;
  double next_step;
  link_data_t *link = &links[link_idx];
  double vol = link_volumes[link_idx];
    
  flowstepsize = link->capacity / (double)NUM_STEPS;
  if (vol <= 0.0) /* should never be -ve but stops warning on == */
    return flowstepsize;
  if (vol >= link->capacity)
    return FLOATINF;
  next_step = ceil(vol / flowstepsize);
  next_step_vol = next_step * flowstepsize;
  if (next_step_vol >= vol)
    return vol;
  else
    return next_step_vol - vol;
}

/* 
 * adjlist_entry_compar() - qsort comparison function for adlist entries
 * 
 * Compares by 'from' node number first  then by 'to' node number if equal
 *
 */
int adjlist_entry_compar(const void *ent1, const void *ent2)
{
  const link_data_t *e1 = (const link_data_t *)ent1;
  const link_data_t *e2 = (const link_data_t *)ent2;
  
  if (e1->init_node < e2->init_node)
    return -1;
  else if(e1->init_node > e2->init_node)
    return 1;
  else
    return ( e1->term_node < e2->term_node ? -1 : 
             (e1->term_node > e2->term_node ? 1 : 0) );
}

/*
 * adjlist_to_packed_arrays() - convert adjlist struct array to packed arrays
 *
 * Packed array format is like
 * like the packed column ("Harwell-Boeing") format used for sparse
 * matrices in some FORTRAN linear algebra routines.
 * For each node i, Va[i] is the first and Va[i+1]-1 the last index into
 * Ea containing the adjacent nodes to i and Wa giving the costs of
 * those edges from node i.
 *
 * Parameters:
 *        adjlist (In/Out) - array of link entries forming adjancey list repr.
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
 * IMPORTANT: The adjlist input must already be sorted.
 */
void adjlist_to_packed_arrays(link_data_t adjlist[],
                                     long num_nodes,
                                     long num_edges,
                                     long Va[],
                                     long Ea[],
                                     double Wa[])
{
  /* sort by 'from' node ascending and within that by 'to' node ascending */
  /* qsort(adjlist, num_edges, sizeof(link_data_t), adjlist_entry_compar); */
  /* no longer need to do that as already sorted after parsing net data */
  int v = 0;  /* index into Va */
  int e = 0;  /* index into Ea and Wa and adjlist */

  for (v = 0; v < num_nodes; v++)
  {
    Va[v] = e;
    while (e < num_edges && adjlist[e].init_node == v)
    {
      Ea[e] = adjlist[e].term_node;
      Wa[e] = adjlist[e].free_flow_time; /* use this as initial cost */
      e++;
    }
    if (v != 0 && Va[v] == e)
      fprintf(stderr, "warning: node %d has no links out\n", v);
  }
  Va[num_nodes] = e;
}


/*
 * update_link_step_costs()
 *
 * update the cost of each link using the cost step function and volume
 *
 * Parameters:
 *    net - net structure parsed from net file
 *    link_volumes - volume on each link
 *    link_costs  (OUT) - cost on each link
 *
 * Return value:
 *    None.
 */
void update_link_step_costs(net_data_t *net, double link_volumes[],
                              double link_costs[])
{
  int k;
  for (k = 0; k < net->num_links; k++)
  {
    link_costs[k] = cost_step_function(&net->links[k], link_volumes[k]);
  }
}

/*
 * demand_sum() - compute total O-D demand
 *
 * Parameters: 
 *    demands -   demands array parsed from trips file
 *                for each origin demands[origin] is an array of demands structs
 *                terminated by one with 0  dest (origins and dests start at 1)
 *    num_zones - number of zones (origins in demands array)
 *
 * Return value:
 *    sum of all O-D demand values
 */
double demand_sum(demand_data_t * const demands[], long num_zones)
{
  long orig,dest,i;
  double total_demand = 0;

  for (orig = 1; orig < num_zones; orig++)
    for (i = 0; (dest = demands[orig][i].dest) != 0; i++)
      total_demand += demands[orig][i].demand;
  return total_demand;
}


/*
 *  def total_link_cost() - 
 *
 *    Compute the total of link costs (volume * costfunction(volume))
 *    over all links at current volumes on those links.
 *    This "measure of effectiveness" is often known as 
 *    "VHT" (Vehicle Hours Travelled).
 *
 *    Parameters:
 *       net - net struct containintg links data
 *       link_volumes - volume on each link
 *       link_costs - cost on each link
 *
 *       
 *
 *    Return value:
 *       total link cost
 */
   
double total_link_cost(const net_data_t *net,
                              const double link_volumes[],
                              const double link_costs[])
{
  int k;
  double total_cost = 0;
  for (k = 0; k < net->num_links; k++)
    total_cost += link_volumes[k]  * link_costs[k];
  return total_cost;
}

/*
 * print_flow() - output flow data giving volume and cost computed by tap
 * 
 * Parameters:
 *     fp - open (write) file pointer to write output to
 *     net - net structure parsed from net file
 *     link_volumes - volume on each link
 *     link_costs  (OUT) - cost on each link
 *
 * Return value:
 *      None.
 */
void print_flow(FILE *fp, const net_data_t *net, 
                       const double link_volumes[],
                       const double link_costs[])
{
  int k;

  fprintf(fp, "<NUMBER OF NODES>\t%ld\n", net->num_nodes);
  fprintf(fp, "<NUMBER OF LINKS>\t%ld\n", net->num_links);
  fprintf(fp, "<END OF METADATA>\n");
  fprintf(fp, "\n");
  fprintf(fp, "~\tTail\tHead\t:\tVolume\tCost\t;\n");
  for (k = 0; k < net->num_links; k++)
    fprintf(fp, "\t%ld\t%ld\t:\t%.15f\t%.15f;\n", net->links[k].init_node,
            net->links[k].term_node, link_volumes[k], link_costs[k]);
}

/*
 * line_search() - line search to find optimal step size by bisection method
 *
 * The 0 <= labmda <=1 that solves 
 *  \Sigma_a \integral_0^{x_a^n} + \alpha (\y_a^n - x_a^n) t_a(\omega) d\omega
 * is found by bisection method (Bolzano method) where
 * the y_a are the (auxilliary) flows (link delta volumes) and x_a the
 * link volumes and t(\omega) is the travel time function (the BPR cost function)
 *
 * Code derived from the LinkSDLineSearch() function in ls_bisect.cpp
 * from  Bar-Gera's sample Frank-Wolfe implemetatnion at
 * http://www.bgu.ac.il/~bargera/tntp/
 *
 * Parameters:
 *       net - net structure parsed from net file
 *      link_volumes - volume on each link
 *      link_delta_volumes - delta volume on each link (search direction)
 *      num_links - number of links (size of link volume vectors)
 *
 * Return value:
 *      step size lambda
 */
double line_search(const net_data_t *net,
                   double link_volumes[], double link_delta_volumes[])
{
  const int MAX_NO_BISECTITERATION = 1000; /* avoid infinite loops */

  int n;
  double lambda_left, lambda_right, lambda;
  double gradient;
  
  gradient=link_directional_derivative(net, link_volumes, link_delta_volumes, 0);
  if (gradient >= 0)
    return 0;
  gradient=link_directional_derivative(net, link_volumes, link_delta_volumes, 1);
  if (gradient <= 0)
    return 1;
  lambda_left = 0;
  lambda_right = 1;
  lambda = 0.5;

  gradient=link_directional_derivative(net,
                                       link_volumes,link_delta_volumes,lambda);
  if (gradient <= 0)
    lambda_left = lambda;
  else
    lambda_right = lambda;
  lambda = 0.5 * (lambda_left + lambda_right);

  for (n = 0; lambda_left == 0 && n < MAX_NO_BISECTITERATION; n++)
  {
    gradient=link_directional_derivative(net,
                                         link_volumes,link_delta_volumes,lambda);
    if (gradient <= 0)
      lambda_left = lambda;
    else
      lambda_right = lambda;
    lambda = 0.5 * (lambda_left + lambda_right);
  }
  return lambda_left;
}


/*
 * apply_trip_modification() - apply a single O-D daemdn data modifcation
 *
 * Parameters:
 *     demands (in/out) - O-D demand data to mofify
 *     mod - trip modificatino to apply to demands
 *     num_zones - number of origins in demands array and dests in each
 *
 * Retrun value:
 *    0 if ok else -1
 */
int apply_trip_modification(demand_data_t *demands[], trip_mod_t *mod,
                             long num_zones)
{
  long i,j;

  /* TODO shoudl sort demand data to use binary search for efficiency instead
     of linear searches */

  if (mod->origin > num_zones || mod->dest > num_zones ||
      mod->origin < 0 || mod->dest < 0)
  {
    fprintf(stderr, "ERROR: bad origin %ld or dest %ld in trip modification\n",
            mod->origin, mod->dest);
    return -1;
    
  }

  switch (mod->trip_modtype)
  {
    case TRIP_MOD_TYPE_POINT_E:
      fprintf(stderr, "[tripchange %s] %ld -> %ld multiplier %f\n",
              mod->trip_change_id, mod->origin, mod->dest, mod->multiplier);
      for (i = 0; demands[mod->origin][i].dest != 0; i++)
      {
        if (demands[mod->origin][i].dest == mod->dest)
        {
          demands[mod->origin][i].demand *= mod->multiplier;
          break;
        }
      }
      break;

    case TRIP_MOD_TYPE_UNIFORM_E:
      if (mod->dest == 0 && mod->origin > 0)
      {
        fprintf(stderr, "[tripchange %s] origin %ld multiplier %f\n",
                mod->trip_change_id, mod->origin, mod->multiplier);
        for (i = 0; demands[mod->origin][i].dest != 0; i++)
          demands[mod->origin][i].demand *= mod->multiplier;
      }
      else if (mod->origin == 0 && mod->dest > 0)
      {
        fprintf(stderr, "[tripchange %s] dest %ld multiplier %f\n",
                mod->trip_change_id, mod->dest, mod->multiplier);
        for (i = 1; i <= num_zones; i++)
        {
          for (j = 0; demands[i][j].dest != 0; j++)
          {
            if (demands[i][j].dest == mod->dest)
              demands[i][j].demand *= mod->multiplier;
          }
        }
      }
      else
      {
        fprintf(stderr, 
                "ERROR: both origin and dest nonzero in uniform trip mod\n");
        return -1;
      }
      break;
      
    case TRIP_MOD_TYPE_ALL_E:
      fprintf(stderr, "[tripchange %s] ALL multiplier %f\n",
              mod->trip_change_id,  mod->multiplier);
      for (i = 1; i <= num_zones; i++)
      {
        for (j = 0; demands[i][j].dest != 0; j++)
          demands[i][j].demand *= mod->multiplier;
      }
      break;
      
    default:
      fprintf(stderr, "bad trip modification type %d\n", mod->trip_modtype);
      return -1;
      break;
  }
  return 0;
}
