#!/usr/bin/env python
###############################################################################
#
# compare_pairwise_linkflow_summary.py - compare pair/sum VHT and link flows
#
# File:    compare_pairwise_linkflow_summary.py
# Author:  Alex Stivala
# Created: May 2011
#
# $Id: compare_pairwise_linkflow_summary.py 628 2011-08-30 00:02:53Z astivala $
# 
###############################################################################

"""
Given the output of the summary.sh script which shows pairwise DeltaVHT
and Sum DeltaVHT values for each pair of changeids, and the actual
link flow outputs for these individual changes, produce sumamries
to compare relative different in VHT between pairwise and sum of individual,
and the overlap in signficiant flow changes on individual links between
the two individual changes, in order to measure correlations
between them (to see if we can use some measure of overlap between 
links with signficnat link flow changes in idnvidual changes to determine
which we should actually have to measure pairwise rather than summing
the indivudal VHT changes).

See usage information in comment documentatino for main.
"""

import sys,os
from time import strftime, localtime
from math import sqrt

from parsetapfiles import parse_net_file,parse_flow_file,parse_node_file,parse_mod_file
from comparelinkflowchanges import compare_delta_flows
from dijkstra import Dijkstra,shortestPath
from traffic_assign import net_to_graph


#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

DELTA_VHT_THRESHOLD = 1000 # ignore abs(DeltaVHT ) < this threshold

#-----------------------------------------------------------------------------
#
# Classes
#
#-----------------------------------------------------------------------------
class Row:
    """
    A dummy class used as struct to represent each row in table
    """
    pass
#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " netfilename nodefilename mods_filename orig_flows_file flow_file_prefix\n")
    sys.exit(1)


def main():
    """
    main for compare_pairwise_linkflow_summary.py

    Usage: compare_pairwise_linkflow_summary.py netfilename nodefilename mods_file orig_file_file flow_file_prefix

      netfilename     is th ename of the network data file defining the road
                       network
      mods_file          is the file containg description of road upgrades
                        (also used as input in tap_frankwolfe_mpi program)
      nodefilename   is name of node file giving x,y cooridinates for nodes
      orig_file_filw is the filename of the file with original (no changes)
                        flow output from tap_frankwfolfe_mpi
      flow_file_prefix is the prefix (path+fileprefix) of flow (output
                       of tap_frankwolfe_mpi) files, to have changeids
                       appended to get link flows for each change.

    Input on stdin is the output of summarize.sh from tap/
       
    Output on stdout is (one line for eaach link with 
     significant changed flow):

     changeid1	changeid2	RelDiff	Overlap NumAugmenting SumFracChagnes	CentroidDistance	normOverlap	normSumFracChanges	NetworkDistance

    Reldiff              - relative difference between pairwise and sum VHT
    Overlap              - number of links with signficiant flow changes
                           (according to compare_delta_flows())
    NumAugmenting        - number of links with signifcant changes in same
                           direction (i.e. both increase or both decrease)
    SumFracChanges       - sum of all the flow chnages in the signficnatly
                           changed links
    CentroidDistance     - Euclidean distance between centroids of upgrade
                           node co-ordinates
    normOverlap          - Overlap normalized to [0,1] by dividing by max
    normSumFracChanges   - SumFracChanges normalized to [0,1] by div by abs max
    NetworkDistance      - Graph distance between upgrade nodes by shortest
                           path using unit edge weights on road network
                           (arbitrarily use lowest node number in set of 
                           upgrade edges as the nodes for this computation)
    
    Reldiff is from summarize.sh and others are compuated from link flow
    files by this script.

    Example usage:

    summarize.sh | compare_pairwise_linkflow_summary.py ChicagoRegional_net.txt ChicagoRegional_node.txt ChicagoRegional_mods.txt ChicagoRegional_flows.txt ChicagoRegional_flows_
    
    """
    if len(sys.argv) != 6:
        usage(os.path.basename(sys.argv[0]))

    netfilename = sys.argv[1]
    nodefilename = sys.argv[2]
    mods_file = sys.argv[3]
    orig_flows_filename = sys.argv[4]
    flowfile_prefix = sys.argv[5]

    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())

    net = parse_net_file(netfilename)
    nodexydict = parse_node_file(nodefilename)
    netmod = parse_mod_file(mods_file)
    orig_flows = parse_flow_file(open(orig_flows_filename))

    for link in net.links:
        link.cost = 1 # unit weights for netgraph
    netgraph = net_to_graph(net)
    
    sys.stdout.write("# Generated by: " + " ".join(sys.argv) + "\n")
    sys.stdout.write("# On: " + timestamp + "\n")
    sys.stdout.write("# From input:\n")


    rows = []
    infh = sys.stdin
    for line in infh:
        if line[0] == "#":
            sys.stdout.write('# ' + line)
            continue
        if line[:9] == "ChangeId1":
            sys.stdout.write("ChangeId1\tChangeId2\tRelDiff\tOverlap\tNumAugmenting\tSumFracChagnes\tCentroidDistance\tnormOverlap\tnormSumFracChanges\tNetworkDistance\tchange1CentroidX\tchange1CentroidY\tchange2CentroidX\tChange2CentroidY\n")
            continue # skip header line
        (changeid1,changeid2,PairDeltaVHT,SumDeltaVHT,RelDiff) = line.split()
        if abs(int(PairDeltaVHT)) < DELTA_VHT_THRESHOLD:
            sys.stderr.write("skipped [%s %s] as %d < %d\n" % 
                              (changeid1,changeid2, abs(int(PairDeltaVHT)),
                               DELTA_VHT_THRESHOLD))
            continue
        change1_flowfilename = flowfile_prefix + changeid1 + ".txt"
        change2_flowfilename = flowfile_prefix + changeid2 + ".txt"
        change1_flows = parse_flow_file(open(change1_flowfilename))
        change2_flows = parse_flow_file(open(change2_flowfilename))
        linkflow_changes = compare_delta_flows(net, orig_flows,
                                               change1_flows, change2_flows)
        change1mods = [mod for mod in netmod if mod.change_id == changeid1]
        change2mods = [mod for mod in netmod if mod.change_id == changeid2]
        change1nodes = ([mod.mod_link.init_node for mod in change1mods] + 
                        [mod.mod_link.term_node for mod in change1mods])
        change2nodes = ([mod.mod_link.init_node for mod in change2mods] + 
                        [mod.mod_link.term_node for mod in change2mods])
        change1centroidX =(sum([nodexydict[node][0] for node in change1nodes])/
                          float(len(change1nodes)))
        change1centroidY =(sum([nodexydict[node][1] for node in change1nodes])/
                          float(len(change1nodes)))
        change2centroidX =(sum([nodexydict[node][0] for node in change2nodes])/
                          float(len(change2nodes)))
        change2centroidY =(sum([nodexydict[node][1] for node in change2nodes])/
                          float(len(change2nodes)))
        centroidDistance = sqrt( (change2centroidX - change1centroidX)**2  +
                                 (change2centroidY - change1centroidY)**2 )

        (pathlen, shortestpath) = Dijkstra(netgraph, min(change1nodes),
                                           min(change2nodes))
        
        if len(linkflow_changes) > 0:
            sum_fracchanges = sum([fracchange1 + fracchange2 for
                                   (initnode,termnode,fracchange1,fracchange2)
                                   in linkflow_changes])
            # number of fractional changes that have the save sign
            # (ie both increase or both decrease flow, rathern than
            #  having opposite effects)
            num_augmenting = len([initnode for
                                   (initnode,termnode,fracchange1,fracchange2)
                                   in linkflow_changes
                                   if cmp(fracchange1,0)==cmp(fracchange2,0)])
        else:
            sum_fracchanges = 0;
            num_augmenting = 0;

        row = Row()
        row.changeid1 = changeid1
        row.changeid2 = changeid2
        row.RelDiff = float(RelDiff)
        row.Overlap = len(linkflow_changes)
        row.NumAugmenting = num_augmenting
        row.SumFracChanges = sum_fracchanges
        row.CentroidDistance = centroidDistance
        row.NetworkDistance = pathlen[min(change2nodes)]
        row.change1centroidX = change1centroidX
        row.change1centroidY =change1centroidY
        row.change2centroidX = change2centroidX
        row.change2centroidY = change2centroidY

        rows.append(row)

    maxOverlap = max([abs(row.Overlap) for row in rows])
    maxSumFracChanges = max([abs(row.SumFracChanges) for row in rows])
    
    for row in rows:    
        try:
            normOverlap = abs(row.Overlap) / float(maxOverlap)
        except ZeroDivisionError:
            normOverlap = float("NaN")
        try:
            normSumFracChanges = abs(row.SumFracChanges) / maxSumFracChanges
        except ZeroDivisionError:
            normSumFracChanges = float("NaN")

        sys.stdout.write("%s\t%s\t% f\t%d\t%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\n" % (
                                            row.changeid1,
                                            row.changeid2,
                                            row.RelDiff,
                                            row.Overlap,
                                            row.NumAugmenting,
                                            row.SumFracChanges,
                                            row.CentroidDistance,
                                            normOverlap,
                                            normSumFracChanges,
                                            row.NetworkDistance,
                                            row.change1centroidX,
                                            row.change1centroidY,
                                            row.change2centroidX,
                                            row.change2centroidY
                                            ))

if __name__ == "__main__":
    main()
    
