#!/usr/bin/env python
###############################################################################
#
# netnodedeltaflow2svg.py - convert traffic model net, node and two sets
#                           of flow data to SVG show flow differences
#
# File:    nodenetdeltaflow2svg.py
# Author:  Alex Stivala
# Created: March 2011
#
# $Id: nodenetdeltaflow2svg.py 273 2011-05-05 06:42:15Z astivala $
# 
###############################################################################

"""
Convert the test data _node, _net, (traffic assignment input)
and two sets of flow (traffic assignment output) data in the format from

http://www.bgu.ac.il/~bargera/tntp/

provided by Hillel Bar-Gera

to SVG format for visualization, showing the differences between the
two sets of flow data.
"""

import sys,os
from time import strftime,localtime
from math import ceil

from parsetapfiles import *

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def write_map_flow_delta_svg(fh, net, xydict, flows1, flows2):
    """
    Write SVG to draw map with differences in flows to open filehandle

    Parameters:
        fh  - open filehandle to write SVG to
        net - class Net object defining nodes and links
        xydict - dict {nodeid:(x,y)} giving positions of nodes
        flows1 - LinkFlow list giving volume and cost on each link for
                 original assignment
        flows2 - LinkFlow list giving volume and cost on each link for
                 new assignment to compare to original

    Return value;
       None.
    """
    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())


    # get link information from net and flows as dict keyed by (from,to)
    flow1_dict = flow_to_linkdict(flows1)
    flow2_dict = flow_to_linkdict(flows2)

    # compute deltas between flows (volumes)
    deltaflow_dict = {}
    for ((nfrom, nto), linkflow) in flow2_dict.iteritems():
        assert(nfrom == linkflow.init_node)
        assert(nto == linkflow.term_node)
        if flow1_dict.has_key((nfrom,nto)): # only if link was also in flows1
            deltaflow_dict[(nfrom, nto)] = (linkflow.volume -
                                            flow1_dict[(nfrom, nto)].volume)
    
    xmin = int(ceil(min([x for (x,y) in xydict.itervalues()])))
    ymin = int(ceil(min([y for (x,y) in xydict.itervalues()])))
    xmax = int(ceil(max([x for (x,y) in xydict.itervalues()])))
    ymax = int(ceil(max([y for (x,y) in xydict.itervalues()])))

    radius = max(xmax-xmin,ymax-ymin) / 1000
    xmin -= radius
    ymin -= radius
    xmax += radius
    ymax += radius
    fh.write('<?xml version="1.0"?>\n')
    fh.write('<svg xmlns="http://www.w3.org/2000/svg" '
     'xmlns:traffic="http://www.csse.unimelb.edu.au/~astivala/traffic.dtd"  '
             'viewBox="%d %d %d %d" >\n' % (xmin, ymin, xmax, ymax))
    fh.write('<traffic:identification xmlns:traffic="http://www.csse.unimelb.edu.au/~astivala/traffic.dtd" creationTime="%s" commandLine="%s" />\n'
             % (timestamp, ' '.join(sys.argv)))


    # plot links
    maxcapacity = max([link.capacity for link in net.links])

    for link in net.links:
        if (link.init_node < net.first_thru_node or
            link.term_node < net.first_thru_node):
            continue  # exclude links to zones (not real links)
        start_xy = xydict[link.init_node]
        end_xy = xydict[link.term_node]
        if deltaflow_dict.has_key((link.init_node,link.term_node)):
            delta_volume = deltaflow_dict[(link.init_node,link.term_node)]
            if flow1_dict[(link.init_node,link.term_node)].volume == 0:
                fracchange = 0 # avoid div by zero
            else:
                fracchange = abs(delta_volume /
                                 flow1_dict[(link.init_node,link.term_node)].volume)
  
  
            # if abs(delta_volume) < 20:
            #     linkcolor = "black"  # no change within threshold
            # elif delta_volume < 0:
            #     linkcolor = "red"  # flow has decreased above threshold
            # else:
            #     linkcolor = "green" # flow has increased above threshold
                
            if fracchange < 0.2 or abs(delta_volume) < 1: # ignore small abs change
                linkcolor = "black"  # no change within threshold
            elif delta_volume < 0:
                linkcolor = "green"  # flow has decreased above threshold
            else:
                linkcolor = "red"    # flow has increased above threshold
            
        else: # it is a new link
            delta_volume = flow2_dict[(link.init_node, link.term_node)].volume
            linkcolor = "purple" # TODO maybe something else to color volratio

        fh.write('    <line traffic:init_node="%d" traffic:term_node="%d" '
                 'traffic:capacity="%f" traffic:delta_volume="%f" '
                 'traffic:volume="%f" '
                 'x1="%f" y1="%f" x2="%f" y2="%f" ' 
                 'stroke-width="%f" '
                 'stroke="%s" />\n'
                 % (link.init_node, link.term_node,
                    link.capacity, delta_volume,
                    flow2_dict[(link.init_node, link.term_node)].volume,
                     start_xy[0],
                    start_xy[1], end_xy[0], end_xy[1], 
                    # draw line width proportional to capacity
                    (link.capacity / maxcapacity) * radius, 
                    linkcolor))
 #   fh.write('  </g>\n')


    fh.write('</svg>\n')


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " netfilename nodefilename orig_flowfilename new_flowfilename\n")
    sys.exit(1)


def main():
    """
    main for nodenetdeltaflow2svg.py

    Usage: nodenetdeltaflow2svg.py netfilename nodefilename orig+flowfilename new_flowfilename

      netfilename    is name of the net file defining node and links
      nodefilename   is name of node file giving x,y cooridinates for nodes
      orig_flowfilename is the flow output file giving volume on links
                          from original assignment
      new_flowfilename is the flow output file giving volume on links from
                            new assignment to compare to original

    Output is SVG on stdout.

    Example usage:
    
    nodenetdeltaflow2svg.py SiouxFalls_net.txt SiouxFalls_node.txt SiouxFalls_flow.txt.orig SiouxFalls_flow.txt.new
    
    """
    if len(sys.argv) != 5:
        usage(os.path.basename(sys.argv[0]))

    netfilename = sys.argv[1]
    nodefilename = sys.argv[2]
    orig_flowfilename = sys.argv[3]
    new_flowfilename = sys.argv [4]

    net = parse_net_file(netfilename)
    nodexydict = parse_node_file(nodefilename)
    flows1 = parse_flow_file(open(orig_flowfilename))
    flows2 = parse_flow_file(open(new_flowfilename))

    write_map_flow_delta_svg(sys.stdout, net, nodexydict, flows1, flows2)


if __name__ == "__main__":
    main()
    
