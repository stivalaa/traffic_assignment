#!/usr/bin/env python
###############################################################################
#
# expand_flows.py - Expand flows from reduced network back to original links
#
# File:    renumber_nodes.py
# Author:  Alex Stivala
# Created: October 2011
#
# $Id: expandflows.py 818 2011-10-27 03:42:32Z astivala $
#
###############################################################################

"""

Take flow output from TAP run on a network reduced by reducenetwork.py
and the linkmap produced by redudcenetwork.py and convert back to flows
file with all links from the original network. Ie the links that were
composed from multiple ilnks in the original are converted back to original
set of links with the flow from the reduced ('shortcut') link put on each
of the original links.

This mapping is in the format:

    newlink oldink1;oldlink2,...

    NB newlink is followede by TAB charadcter then oldink list
    where newlink and the oldlinks are fromnode,tonode tuples
    (i.e. node numbers delimited by comma) and the list of
    comma-delimited oldlink tuples is delmited by semcolons e.g.

    (11295, 11141)	(11295, 11276);(11276, 11142);(11142, 11141)

Note that these node numbers are all in the OLD (original node
numbering), unlike the output net and node files which have them
renumbered. hence the flow input into this script should not come
directly from the TAP output, but be converted from that back to
original node numbers with renumber_nodes.py



"""


import sys,os
from time import strftime, localtime

from parsetapfiles import parse_flow_file,LinkFlow,flow_to_linkdict,parse_net_file,Link,net_to_linkdict
from traffic_assign import print_flow
from reducenetwork import linkcmp

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def parse_linkmap(fh):
    """
    Read the mapping between new (added) links and the multiple (now
    removed) links they correspdond to is written to
    out_linkmapfilename in the format (one mapping perline):

    newlink oldink1;oldlink2,...

    NB newlink is speatae from oldlink list by a TAB chardcter
    where newlink and the oldlinks are fromnode,tonode tuples
    (i.e. node numbers delimited by comma) and the list of
    comma-delimited oldlink tuples is delmited by semcolons e.g.

    (11295, 11141)	(11295, 11276);(11276, 11142);(11142, 11141)

   Note that these node numbers are all in the OLD (original node
   numbering), unlike the output net and node files which have them
   renumbered.

    Return value: dict { newlink_tuple : [list of oldlink tuples }
       mapping a new link (as (fromnode,tonode tuple)) back to list of
       oldlink tuples.
    """
    linkmap = {}
    for line in fh:
        if line[0] == '~':
            continue
        sline = line.split('\t')
        if len(sline) != 2:
            sys.stderr.write('ERROR: bad line "%s"\n' % (line))
            return None
        newlink = eval(sline[0])
        oldlink_list = [eval(t) for t in sline[1].split(';')]
        if linkmap.has_key(newlink):
            sys.stderr.write("ERROR duplicate link key %d\n" % (newlink))
            return None
        linkmap[newlink] = oldlink_list
    return linkmap


def bpr_cost_function(link, Q):
    """
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
    """
    alpha = link.B
    beta  = link.power
    Q_max = link.capacity
    assert(Q_max > 0);
    return link.free_flow_time * (1.0 + alpha * pow(Q / Q_max, beta))
#  /* return link->free_flow_time * (1.0 + alpha *  */
#  /*                                (link->B / pow(link->capacity, link->power)) *  */
#  /*                                pow(Q, beta)); */


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    Print usage message and exit
    """

    sys.stderr.write("Usage: " +progname + " linkmap_file\n")
    sys.exit(1)


def main():
    """
    main for expandflows.py

    Usage: expandflows.py net_file linkmap_file 

    net_file is the network file for the full (expanded, not reduced) network
    linkmap_file is the reduced network link mapping from from reducenetwork.py

    Input is a flow file on stdin (renumberd from reduced network TAP
    output with rennumber_nodes.py) , output is flow file with
    links expanded back to original, with flows on them from reduced network TAP solution,
    on stdout.
   
    Example usage:

        renumber_nodes.py melbourne_reduce_number_map.txt  <  Melbourne_reduced_flows.txt \
        | expandflows.py melbourne_net.txt melbourne_reduced_linkmap.txt  > Melbourne_expanded_flows.txt

    """
    if len(sys.argv) != 3:
        usage(os.path.basename(sys.argv[0]))

    net_file = sys.argv[1]
    linkmap_file = sys.argv[2]

    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
    net = parse_net_file(net_file)
    linkdict = net_to_linkdict(net)
    linkmap = parse_linkmap(open(linkmap_file))
    
    flows = parse_flow_file(sys.stdin)
    flow_dict = flow_to_linkdict(flows)
    expanded_flow_dict = {}
    for linkflow in flow_dict.itervalues():
        if linkmap.has_key((linkflow.init_node,linkflow.term_node)):
            # add all the old links that this link coaeseced into expanded dict
            for (fromnode, tonode) in linkmap[(linkflow.init_node,linkflow.term_node)]:
                lf = LinkFlow()
                lf.init_node = fromnode
                lf.term_node = tonode
                lf.volume = linkflow.volume
                lf.cost = bpr_cost_function(linkdict[(lf.init_node,lf.term_node)], lf.volume)
                expanded_flow_dict[(lf.init_node, lf.term_node)] = lf
        else:
            # not an expanded link,just copy the original linkflow to
            # expanded linkflow dict
            lf = LinkFlow()
            lf.init_node = linkflow.init_node
            lf.term_node = linkflow.term_node
            lf.volume = linkflow.volume
            lf.cost = linkflow.cost
            expanded_flow_dict[(lf.init_node, lf.term_node)] = lf

    outfh = sys.stdout
    outfh.write('~ Generated by: ' + ' '.join(sys.argv) + '\n')
    outfh.write('~ On: ' + timestamp + '\n')

    # dodgy: fake net object just t oget num_nodes and num_keys used in print_flow
    class Net:
        pass
    net = Net()
    net.num_nodes = max([max(fromnode,tonode) for (fromnode,tonode) in expanded_flow_dict.iterkeys()])
    net.num_links = len(expanded_flow_dict)

    print_flow(outfh, net, expanded_flow_dict)

    
if __name__ == "__main__":
    main()

