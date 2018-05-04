#!/usr/bin/env python
###############################################################################
#
# traffic_assign.py - traffic assignment by greedy method
#
# File:    traffic_assign.py
# Author:  Alex Stivala
# Created: March 2011
#
# OBSOLETE - use C implementation in tap/ now
#
# $Id: traffic_assign.py 797 2011-10-11 00:32:16Z astivala $
# 
###############################################################################

"""
Traffic assignment by greedy algorithm.

"""

import sys,os
from time import strftime,localtime,clock
from math import ceil,floor

from dijkstra import Dijkstra,shortestPath
from parsetapfiles import *


#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# steps in cost step function
NUM_STEPS = 3.0


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def total_link_cost(net):
    """
    Compute the total of link costs (volume * costfunction(volume))
    over all links at current volumes on those links

    Parameters:
       net - Net object as returned by parse_net_file()

    Return value:
       total link cost
    """
    return sum([link.volume * link.cost for link in net.links])

def cost_step(link, volume):
    """
    Given a Link object and volume on that link, return the step
    (in the division of volume up to capacity) that the volume is at
    step function at that volume

    Parameters:
        link - Link object
        volume- volume to get cost for on link

    Return value:
        step (0, 1, .. NUM_STEPS) that the volume is at on that link
    """
    if volume >= link.capacity:
        return  NUM_STEPS
    flowstepsize = link.capacity / NUM_STEPS
    step = floor(volume / flowstepsize)
    return step

def distance_from_next_step(link):
    """
    Given a Link object return the distance (volume difference) from
    the point where it will go up to the next step in the step cost function

    Parameters:
        link - Link object

    Return value:
        distance (additional volume required) to go to next step in cost step fn
    """
    flowstepsize = link.capacity / NUM_STEPS
    if link.volume == 0.0:
        return flowstepsize
    elif link.volume >= link.capacity:
        return float("inf")
    next_step = ceil(link.volume / flowstepsize)
    next_step_vol = next_step * flowstepsize
    if next_step_vol == link.volume:
        return next_step_vol
    else:
        return next_step_vol - link.volume
    
def cost_step_function(link, volume):
    """
    Given a Link object and volume on that link, return the value of cost
    step function at that volume

    Parameters:
        link - Link object
        volume- volume to get cost for on link

    Return value:
        cost at volume on link
    """
    step = cost_step(link, volume)
    coststepsize = 2.0 # FIXME 
    if volume > link.capacity:
        cost = link.free_flow_time  + step * coststepsize # TODO something else?
    else:
        cost = link.free_flow_time  + step * coststepsize

    return cost
    


def assign_all_or_nothing(net, demandsdict, linkdict):
    """
    Assign flow from O-D demands to links according to current shortest
    paths (lowest costs) on netgraph
    (all-or-nothing assignment)

    Parameters:
        net - Net object as returned by parse_net_file()
        demandsict -  O-D demand data represented as {(orig,dest):demand}
        linkdict (in/out) - dict {(from,to):Link} of Link objects with volume
                   attribute assigned according to shortest paths


    Return value:
    """
    netgraph = net_to_graph(net)
    for ((nfrom,nto),link) in linkdict.iteritems():
        link.volume = 0.0
        
    for ((orig, dest), routeflow) in demandsdict.iteritems():
        pathnodes = shortestPath(netgraph, orig, dest)
        path = zip(pathnodes, pathnodes[1:]) #convert to list of edge tuples
        for (nfrom, nto) in path:
            linkdict[(nfrom, nto)].volume += routeflow
            

def get_shortest_path(pred, orig, dest):
    """
    Get the shortest path to dest given the predecessor list pred
    from Dijkstra's algorithm

    Parameters:
       pred - list of predecessor nodes from Dijkstra's algorithm with orig 
              as source
       source - source node
       dest  - destinatino node

    Return value:
      list of nodes on shortest path to dest from orig
    """
    path = []
    v = dest
    while v != orig:
        path.append(v)
        v = pred[v]
    path.append(orig)
    path.reverse()
    return path


def assign(net, demands, demandsdict, linkdict):
    """
    Assign flow from O-D demands to links according to current shortest
    paths (lowest costs) on netgraph
    Assigns flow to path only up to the amount of flow that
    takes the first link to do so up to the next step in cost function,
    subtracting that assigned flow from the O-D demand using that path

    Parameters:
        net - Net object as returned by parse_net_file()
        demands (in/out) - dict of dicts { origin : { destination : demand } }
                  from parse_trips_file()
        demandsict (in/out) - O-D demand data represented as {(orig,dest):demand}
                   from demands_to_demands_dict(demands)
        linkdict (in/out) - dict {(from,to):Link} of Link objects with volume
                   attribute assigned according to shortest paths


    Return value:
    """
    netgraph = net_to_graph(net)
    for (orig, demand_dict) in demands.iteritems():
#        sys.stderr.write("running dijjkstra for orig %d..." % orig)
        t0 = clock()
        (d, pred) = Dijkstra(netgraph, orig)
        t1 = clock()
#        sys.stderr.write("done (%f ms).\n" % ((t1 - t0)*1000))
        for (dest, routeflow) in demand_dict.iteritems():
            if orig == dest or routeflow == 0.0:
                continue
            pathnodes = get_shortest_path(pred, orig, dest)
            path = zip(pathnodes, pathnodes[1:]) #convert to list of edge tuples
            path_links = [linkdict[from_to_tuple] for from_to_tuple in path]
            # assign flow to path only up to the amount of flow that
            # takes the first link to do so up to the next step in cost function
            addflow = min(routeflow,
                min([distance_from_next_step(link) for link in path_links]))
            for link in path_links:
                link.volume += addflow
            demandsdict[(orig, dest)] -= addflow
            demand_dict[dest] -= addflow 


def update_link_costs(net):
    """
    Update the cost element in each Link with using cost step function
    and volume
    """
    for link in net.links:
        link.cost = cost_step_function(link, link.volume)


def print_flow(fh, net, linkdict):
    """
    Output the flow data giving volume and cost from output of traffic assign
    
    Paraemters:
        fh - open (write) filehande to write flow data to
        net - Net object as returned by parse_net_file()
        linkdict - dict {(from,to):Link} of Link objects with volume
                   attribute assigned according to shortest paths

    Return value:
       None.

       Output  volume and cost on each link in format like:

       <NUMBER OF NODES>       25
       <NUMBER OF LINKS>       76
       <END OF METADATA>

       ~       Tail    Head    :       Volume  Cost    ;
               1       2       :       4494.5008499645437041   6.0008161234622576785   ;
               1       3       :       8119.1900669362912595   4.0086912217072878661  

    """
    assert(net.num_links == len(list(linkdict.iterkeys())))
    fh.write("<NUMBER OF NODES>\t%d\n" % net.num_nodes)
    fh.write("<NUMBER OF LINKS>\t%d\n" % net.num_links)
    fh.write("<END OF METADATA>\n")
    fh.write("\n\n")
    fh.write("~\tTail\tHead\t:\tVolume\tCost\t;\n")

    # sort by fromNode then toNode for output
    linklist = list(linkdict.iteritems())
    linklist.sort()

    for ((nfrom, nto), link) in linklist:
        assert(link.init_node == nfrom)
        assert(link.term_node == nto)
        fh.write("\t%d\t%d\t:\t%.15f\t%.15f\t;\n" %
                 (link.init_node, link.term_node, link.volume, link.cost))
           
    
#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " netfilename  demandfilename\n")
    sys.exit(1)


def main():
    """
    main for traffic_assign.py

    Usage: traffic_assign.py netfilename  demandfilename

      netfilename    is name of the net file defining node and links
      demandfilename is name of Origin-Destination demand file

      Output is link flows on stdout.
      
    Example usage:
    
    traffic_assign.py SiouxFalls_net.txt SiouxFalls_trips.txt
    
    """
    if len(sys.argv) != 3:
        usage(os.path.basename(sys.argv[0]))

    netfilename = sys.argv[1]
    demandfilename = sys.argv[2]

    net = parse_net_file(netfilename)
    demands = parse_trips_file(demandfilename)

    linkdict = net_to_linkdict(net)
    demandsdict = demands_to_demanddict(demands)

    total_od_flow = sum(demandsdict.itervalues())
    sys.stderr.write("total OD flow = %f\n" % total_od_flow)

    voldict = dict([ ((link.init_node, link.term_node), link.volume)
                     for link in net.links])

    iternum = 0
    while max(demandsdict.itervalues()) > 0.0:

        old_voldict = dict(voldict)
        assign(net, demands, demandsdict, linkdict)
        update_link_costs(net)
                
        voldict = dict( [((link.init_node, link.term_node), link.volume)
                         for link in net.links])

        iternum += 1

        delta_link_vol =dict([((link.init_node, link.term_node),
                                  voldict[(link.init_node, link.term_node)]
                                  - old_voldict[(link.init_node, link.term_node)])
                                  for link in net.links])
        delta_link_costs =  [ delta_link_vol[(link.init_node, link.term_node)]
                              * link.volume for link in net.links ]
        
        avg_excess_cost = sum(delta_link_costs) / total_od_flow

        sys.stderr.write("iter = %d total unsatisfied demand = %f avg excess cost = %f\n" % 
                         (iternum, sum(demandsdict.itervalues()),
                          avg_excess_cost))



    
    print_flow(sys.stdout, net, linkdict)


if __name__ == "__main__":
    main()
