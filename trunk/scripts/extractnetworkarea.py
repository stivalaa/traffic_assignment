#!/usr/bin/env python
##############################################################################
#
# extractnetworkarea.py - Extract part of a network
#
# File:    extractnetworkarea.py
# Author:  Alex Stivala
# Created: October 2011
#
# $Id: extractnetworkarea.py 810 2011-10-25 03:39:18Z astivala $
#
#
##############################################################################

"""
Extract part of a network inside a specified rectangular area

Usage:
     extractnetworkarea netfilename nodefilename out_netfilename out_nodefilename x1 y1 x2 y2

      netfilename    is name of the net file defining node and links
      nodefilename   is name of node file giving x,y cooridinates for nodes
      out_netfilename is name of output net file for extraced netwirk
      out_nodefilename is name of output node file for extracted network
      x1 y1 co-ordinates of upper left corner of rectangle to extract
      x2 y2 co-ordinates of lower right corner of rectangle to extract

    Since node numbers have to be sequential (from 1), a mapping from
    original node numbers to renumbered (after deleting nodes outside
    rectangle) sequential nuode numbers is written to stdout in form
   
    oldnum newnum

    one per line.
   
"""

import sys
import getopt
from time import strftime, localtime
from parsetapfiles import Link,Net,parse_net_file,parse_node_file,net_to_linkdict

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname +
                     " netfilename nodefilename out_netfilename out_nodefilename x1 y1 x2 y2\n")
    sys.exit(1)

def main():
    """
    See usage message in module header block
    """
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        usage(sys.argv[0])

    if len(args) != 8:
        usage(sys.argv[0])


    net_filename = args[0]
    nodefilename = args[1]
    out_net_filename = args[2]
    out_node_filename = args[3]
    x1 = float(args[4])
    y1 = float(args[5])
    x2 = float(args[6])
    y2 = float(args[7])

    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())

    net = parse_net_file(net_filename)
    nodexydict = parse_node_file(nodefilename)

    num_nodes = net.num_nodes

#    sys.stderr.write(str([ (nodexydict[node][0], nodexydict[node][1]  ,
#                       nodexydict[node][0], nodexydict[node][1])  for node  in range(3000,3010)]))
    
    extracted_nodes = [node for node in xrange(1, num_nodes+1) if
                       nodexydict[node][0] >= x1 and nodexydict[node][1] >= y1 and
                       nodexydict[node][0] <= x2 and nodexydict[node][1] <= y2]

    num_extracted_nodes = len(extracted_nodes)

    extracted_nodes_dict = dict([(node, True) for node in extracted_nodes])
    extracted_links = [link for link in net.links if
                       extracted_nodes_dict.has_key(link.init_node) and
                       extracted_nodes_dict.has_key(link.term_node)]

    sys.stderr.write("extracted %d nodes and %d links inside rectangle\n" %
                     (num_extracted_nodes, len(extracted_links)))

    removed_nodes = [node for node in xrange(1, num_nodes+1) if
                     node not in extracted_nodes_dict]

    # now remove nodes that no longer have any links in or out.
    # first create dict { node : True } for nodes that have a link out AND
    # one for links in
    node_linksout = {}
    for link in extracted_links:
        node_linksout[link.init_node] = True
    node_linksin = {}
    for link in extracted_links:
        node_linksin[link.term_node] = True

    nodes_with_no_links = [node for node in extracted_nodes if not node_linksout.has_key(node)
                           and not node_linksin.has_key(node)]
#    sys.stderr.write('xxxx ' + str(nodes_with_no_links))
    removed_nodes += nodes_with_no_links
    extracted_nodes = list(set(extracted_nodes) - set(nodes_with_no_links))
    num_extracted_nodes = len(extracted_nodes)
    extract_num_zones = len([node for node in extracted_nodes if node < net.first_thru_node])
    
    # But node numbering has to be sequential so build a dict
    # { node_num : new_node_num } mapping the current node numbers,
    # for nodes with links only, to a new sequential numbering
    removed_nodes_dict = dict([(node, True) for node in removed_nodes])
    seqnum = 1
    renumber_dict = {}
    for node in xrange(1, net.num_nodes+1):
        if not removed_nodes_dict.has_key(node):
            renumber_dict[node] = seqnum
            seqnum += 1

    reverse_renumber_dict = dict([(newnum, oldnum) for
                                  (oldnum, newnum) in
                                  renumber_dict.iteritems()])
                                 
    outfh = open(out_net_filename, "w")
    outfh.write('~ Generated by: ' + ' '.join(sys.argv) + '\n')
    outfh.write('~ On: ' + timestamp + '\n')
    outfh.write("<NUMBER OF ZONES> %d\n" % extract_num_zones)
    outfh.write("<FIRST THRU NODE> %d\n" % (extract_num_zones + 1))
    outfh.write("<NUMBER OF NODES> %d\n" % num_extracted_nodes)
    outfh.write("<NUMBER OF LINKS> %d\n" % len(extracted_links))
    outfh.write("<END OF METADATA>\n")
    outfh.write('~ from \tto\tcapacity\tlength\tftime\tB\tpower\tspeed\ttoll\ttype\n')

    for link in extracted_links:
        outfh.write("%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t;\n" %
                    (renumber_dict[link.init_node],
                     renumber_dict[link.term_node],
                     link.capacity,
                     link.length,link.free_flow_time,
                     link.B,link.power,link.speed_limit,
                     link.toll,link.linktype))

    outfh.close()

    node_outfh = open(out_node_filename, "w")
    node_outfh.write("Node\tX\tY\t;\n")
    for node in xrange(1,  num_extracted_nodes+1):
        (x, y) = nodexydict[reverse_renumber_dict[node]]
        node_outfh.write("%d\t%f\t%f\t;\n" % (node, x, y))
    node_outfh.close()

    # write the map file mapping orig numbers to new sequential node numbers
    mapfile_fh = sys.stdout
    mapfile_fh.write('~ Generated by: ' + ' '.join(sys.argv) + '\n')
    mapfile_fh.write('~ On: ' + timestamp + '\n')
    mapfile_fh.write('~oldNodeNum\tnewNodeNum\n')
    for (oldnum,newnum) in sorted(renumber_dict.iteritems()):
        mapfile_fh.write('%d\t%d\n' % (oldnum, newnum))

    
if __name__ == "__main__":
    main()
