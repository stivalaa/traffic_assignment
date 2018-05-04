#!/usr/bin/env python
###############################################################################
#
# flow2rmat.py - convert flow output file to rmat.txt format for SSSP tests
#
# File:    flow2rmat.py
# Author:  Alex Stivala
# Created: September 2011
#
# $Id: flow2rmat.py 825 2011-10-31 06:27:16Z astivala $
# 
###############################################################################

"""
Conver the flow output file from TAP in TAP
format (http://www.bgu.ac.il/~bargera/tntp/) to rmat.txt matrix format
where each line give from node, to node, and cost,
for use the the cuda_sssp/ test programs.

Reads flow file on stdin, write to stdout.
"""

import sys,os

from parsetapfiles import parse_flow_file,LinkFlow


if len(sys.argv) > 1:
    sys.stderr.write("usage: flow2rmat.py < flowfile > rmat.txt.output\n")
    sys.exit(1)

link_flows = parse_flow_file(sys.stdin)
num_nodes = max([max(link.init_node, link.term_node) for link in link_flows])
sys.stdout.write("%d\t%d\t%d\n" % (num_nodes, num_nodes, len(link_flows)))
for link in link_flows:
    sys.stdout.write("%d\t%d\t%.15f\n" % (link.init_node, link.term_node, link.cost))


