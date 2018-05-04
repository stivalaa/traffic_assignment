#!/usr/bin/env python
###############################################################################
#
# comparelinkflowchanges.py - show links with delta flows in common in 2 mods
#
# File:    comparelinkflowchanges.py
# Author:  Alex Stivala
# Created: May 2011
#
# $Id: comparelinkflowchanges.py 352 2011-06-08 05:17:32Z astivala $
# 
###############################################################################

"""
Given two flow output files from two different network changes and the
original flow output file from the unmodified network, show links
that have signficant changes in flow in both modifications.
"""

import sys,os
from math import ceil

from parsetapfiles import *

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def compare_delta_flows(net, orig_flows, change1_flows, change2_flows):
    """
    Compare both change's flows to original and find links where both
    network changes have made signficiant differences in flow.

    Parameters:
        net - class Net object defining nodes and links
        orig_flows - LinkFlow list giving volume and cost on each link for
                    original assignment
        change1_flows - LinkFlow list giving volume and cost on each link for
                     assignment from a change
        change2_flows - LinkFlow list giving volume and cost on each link for
                     assignment from a second change

    Return value;
       list of (initnode, termnode, fracchange1, fracchange2) tuples
       For a new (added) link, fracchange1 and fracchange2 may be None
    """
    changelist = []

    # get link information from net and flows as dict keyed by (from,to)
    orig_flow_dict = flow_to_linkdict(orig_flows)
    change1_flow_dict = flow_to_linkdict(change1_flows)
    change2_flow_dict = flow_to_linkdict(change2_flows)

    # compute deltas between flows (volumes)
    deltaflow1_dict = {}
    for ((nfrom, nto), linkflow) in change1_flow_dict.iteritems():
        assert(nfrom == linkflow.init_node)
        assert(nto == linkflow.term_node)
        if orig_flow_dict.has_key((nfrom,nto)):
            deltaflow1_dict[(nfrom, nto)] = (linkflow.volume -
                                            orig_flow_dict[(nfrom, nto)].volume)
    deltaflow2_dict = {}
    for ((nfrom, nto), linkflow) in change2_flow_dict.iteritems():
        assert(nfrom == linkflow.init_node)
        assert(nto == linkflow.term_node)
        if orig_flow_dict.has_key((nfrom,nto)):
            deltaflow2_dict[(nfrom, nto)] = (linkflow.volume -
                                            orig_flow_dict[(nfrom, nto)].volume)

    for link in net.links:
        if (link.init_node < net.first_thru_node or
            link.term_node < net.first_thru_node):
            continue  # exclude links to zones (not real links)
        if (orig_flow_dict.has_key((link.init_node,link.term_node))):
            delta1_volume = deltaflow1_dict[(link.init_node,link.term_node)]
            delta2_volume = deltaflow2_dict[(link.init_node,link.term_node)]
            if orig_flow_dict[(link.init_node,link.term_node)].volume == 0:
                fracchange1 = 0 # avoid div by zero
                fracchange2 = 0 # avoid div by zero
            else:
                fracchange1 = (delta1_volume /
                        orig_flow_dict[(link.init_node,link.term_node)].volume)
                fracchange2 = (delta2_volume /
                        orig_flow_dict[(link.init_node,link.term_node)].volume)
            if (not ((abs(fracchange1) < 0.2 or abs(delta1_volume) < 1) or
                     (abs(fracchange2) < 0.2 or abs(delta2_volume) < 1)) ):
                changelist.append((link.init_node, link.term_node,
                                   fracchange1, fracchange2))
            else:
                 pass
#                print 'xxx',link.init_node,link.term_node,fracchange1,fracchange2,delta1_volume,delta2_volume
        else: # it is a new link
            # This won't happen anyway when the net file is the net file
            # for the original flows, since the iteration is only over links
            # already  in that net (orig) file
            changelist.append((link.init_node, link.term_node, None, None))

    return changelist



#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " netfilename orig_flowfilename change1_flowfilename change2_flowfilename\n")
    sys.exit(1)


def main():
    """
    main for comparelinkflowchanges.py

    Usage: comparelinkflowchange.py netfilename orig_flowfilename change1_flowfilename change2_flowfilename

      netfilename    is name of the net file defining node and links
      orig_flowfilename is the flow output file giving volume on links
                          from original assignment
      change1_flowfilename is the flow output file giving volume on links from
                            first change
      change2_flowfilename is the flow output file giving volume on links from
                            second change

    Output on stdout is (one line for eaach link with 
     significant changed flow):

     initnode termnode francchange1 fracchange2

    Example usage:

./comparelinkflowchanges.py ~/traffic_assignment/trunk/testdata/ChicagoRegional/ChicagoRegional_net.txt /var/tmp/ChicagoRegional_flows.txt  /var/tmp/ChicagoRegional_flows_03_96_0024.txt  /var/tmp/ChicagoRegional_flows_07_06_0014.txt
    
    
    """
    if len(sys.argv) != 5:
        usage(os.path.basename(sys.argv[0]))

    netfilename = sys.argv[1]
    orig_flowfilename = sys.argv[2]
    change1_flowfilename = sys.argv [3]
    change2_flowfilename = sys.argv [4]

    net = parse_net_file(netfilename)
    orig_flows = parse_flow_file(open(orig_flowfilename))
    change1_flows = parse_flow_file(open(change1_flowfilename))
    change2_flows = parse_flow_file(open(change2_flowfilename))

    changelist = \
    compare_delta_flows(net,orig_flows,change1_flows,change2_flows)
    for (initnode,termnode,fracchange1,fracchange2) in changelist:
        if fracchange1 and fracchange2:
            sys.stdout.write("%d\t%d\t%f\t%f\n" % (initnode, termnode,
                                               fracchange1, fracchange2))
        else:
            sys.stdout.write("%d\t%d\tNEWLINK\n" %(link.init, link.term))

if __name__ == "__main__":
    main()
    
