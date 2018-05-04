#!/usr/bin/env python
###############################################################################
#
# nodenetflow2svg.py - convert traffic model net, node and flow data to SVG
#
# File:    nodenetflow2svg.py
# Author:  Alex Stivala
# Created: March 2011
#
# $Id: nodenetflow2svg.py 812 2011-10-25 23:40:05Z astivala $
# 
###############################################################################

"""
Convert the test data _node, _net, (traffic assignment input)
and flow (traffic assignment output) data in the format from

http://www.bgu.ac.il/~bargera/tntp/

provided by Hillel Bar-Gera

to SVG format for visualization.

"""

import sys,os
import getopt
from time import strftime,localtime
from math import ceil

from parsetapfiles import *

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def write_map_svg(fh, net, xydict, flows, plot_zones = False,
                  plot_centroid_connectors = False,
                  netmod = None,
                  debuglinktype999 = False):
    """
    Write SVG to draw map to open filehandle

    Parameters:
        fh  - open filehandle to write SVG to
        net - class Net object defining nodes and links
        xydict - dict {nodeid:(x,y)} giving positions of nodes
        flows - LinkFlow list giving volume and cost on each link
        plot_zones - if True, plot zone centroids as purple circles
        plot_centroid_connectors - if True, plot centroid connectors as
                                    aqua lines.
        netmod - NetMod object with added/changed roads. If not None,
                 plot these in red.
        debuglinktype999 - if True, colour link type 999 magenta, used for
                 debugging reducenetwork.py

    Return value;
       None.
    """
    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())

    # if upgrades file provided, apply upgardes to road network
    for link in net.links:
        link.change_id = ""
        link.project_cost = ""
    if netmod != None:
        linkdict = net_to_linkdict(net)
        addlinkdict = {} # dict { (initnode,termnod) : True } to check for duplicates
        for mod in netmod:
            if mod.modtype == MODTYPE_ADD:
                if (linkdict.has_key((mod.mod_link.init_node,mod.mod_link.term_node))
                    or addlinkdict.has_key((mod.mod_link.init_node,
                                            mod.mod_link.term_node))):
                    sys.stderr.write('warning: add %s initnode %d termnode %d '
                                     'ignored as duplicate\n' %
                                     (mod.change_id,mod.mod_link.init_node,
                                      mod.mod_link.term_node))
                    continue
                mod.mod_link.change_id = mod.change_id
                mod.mod_link.project_cost = mod.project_cost
                net.links.append(mod.mod_link)
                addlinkdict[(mod.mod_link.init_node,mod.mod_link.term_node)] = True
                sys.stderr.write('[add %s] initnode %d termnode %d capacity %f\n' %
                                 (mod.change_id, mod.mod_link.init_node,
                                  mod.mod_link.term_node, mod.mod_link.capacity))
            elif mod.modtype == MODTYPE_CHANGE:
                for i in xrange(len(net.links)):
                    link = net.links[i]
                    if (link.init_node == mod.mod_link.init_node and
                        link.term_node == mod.mod_link.term_node):
                        sys.stderr.write('[change %s] initnode %d termnode %d '
                                         'old capacity %f new capacity %f\n' %
                                         (mod.change_id, mod.mod_link.init_node,
                                          mod.mod_link.term_node,
                                          link.capacity, mod.mod_link.capacity))
                        net.links.pop(i)
                        mod.mod_link.change_id = mod.change_id
                        mod.mod_link.project_cost = mod.project_cost
                        net.links.append(mod.mod_link)
                        break
            else:
                raise ValueError('unknown modtype ' + str(mod.modtype))


    # get link information from net and flows as dict keyed by (from,to)
#(not needed -iterating over links)    link_dict = net_to_linkdict(net)
    flow_dict = None
    if flows:
        flow_dict = flow_to_linkdict(flows)
    
    
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


    # plot nodes
    # for (nodeid, (x, y)) in xydict.iteritems():
    #     fh.write('    <circle traffic:nodeid="%d" cx="%f" cy="%f" r="%f"  />\n' 
    #              % (nodeid, x, y, radius))

    # plot 'zones' (origin or destination nodes)
    if plot_zones:
        for nodeid in xrange(1, net.num_nodes):
            try:
                (x, y) = xydict[nodeid]
            except KeyError:
                sys.stderr.write("warning: no node " + str(nodeid)+"\n")
                continue
            if (nodeid < net.first_thru_node):
                fh.write('    <circle traffic:nodeid="%d" fill="purple" cx="%f" cy="%f" r="%f"  />\n' 
                         % (nodeid, x, y, radius)) # FIXME this makes map almost unreadable but if any smaller inkscape won't let me select them
            else:
                continue # XXX only plot zones
                fh.write('    <circle traffic:nodeid="%d" fill="red" cx="%f" cy="%f" r="%f"  />\n' 
                         % (nodeid, x, y, radius/10))
    

    # plot links
#    linewidth = radius
    maxcapacity = max([link.capacity for link in net.links if
                       link.init_node >=  net.first_thru_node and
                        link.term_node >= net.first_thru_node])
    sys.stderr.write("max non centroid connector capacity is %f\n" % maxcapacity)
    if maxcapacity >= 999999.000000: # hack to work on Berlin data
        maxcapacity = max([link.capacity for link in net.links if
                           link.init_node >=  net.first_thru_node and
                           link.term_node >= net.first_thru_node and
                           link.capacity < 999999.000000]) 
        sys.stderr.write("ignoring capacity >= 999999.000000, max capacity is now %f\n" % maxcapacity)
#    fh.write('  <g stroke-width="%d">\n' % (linewidth))
    for link in net.links:
        start_xy = xydict[link.init_node]
        end_xy = xydict[link.term_node]

        if  ( (link.init_node < net.first_thru_node or 
                  link.term_node < net.first_thru_node) ):
            if plot_centroid_connectors:
                linkcolor = "aqua"
            else:
                continue
            if flow_dict:
                volume = flow_dict[(link.init_node,link.term_node)].volume
            else:
                volume = 0
        else:
            if flow_dict:
                volume = flow_dict[(link.init_node,link.term_node)].volume
                if link.capacity == 0:
                    volratio = 1
                else:
                    volratio = volume / link.capacity
                if volratio > 1.00:
                    linkcolor = "red"
                elif volratio > 0.95:
                    linkcolor = "yellow"
                elif volratio < 0.05:
                    linkcolor = "black"
                else:
                    linkcolor = "green"
            else:
                volume = 0
                linkcolor = "black"
            if netmod:  # highlight changed/added links
                if link.change_id != "": # was not in original network
                    linkcolor = "red"
                

        # draw line width proportional to capacity
        if link.capacity == 0:
            sys.stderr.write("warning: 0 capacity on link from %d to %d\n" %
                             (link.init_node, link.term_node))
            linewidth = 0.1 * radius
        elif link.capacity >= 999999:
            sys.stderr.write("warning: 999999 capacity on link from %d to %d\n" %
                             (link.init_node, link.term_node))
            linewidth = radius
        else:
            linewidth = (link.capacity / maxcapacity) * radius


        # colour composed links with linktype 999 (used for debugging reducenetwork.py)
        if debuglinktype999 and link.linktype == 999:
            linkcolor = "magenta"

        fh.write('    <line traffic:init_node="%d" traffic:term_node="%d" '
                 'traffic:capacity="%f" traffic:volume="%f" '
                 'traffic:linktype="%d" '
                 'traffic:free_flow_time="%f" '
                 'traffic:B="%f" '
                 'traffic:power="%f" '
                 'traffic:length="%f" '
                 'traffic:speedlimit="%f" '
                 'traffic:change_id="%s" '
                 'x1="%f" y1="%f" x2="%f" y2="%f" ' 
                 'stroke-width="%f" '
                 'stroke="%s" />\n'
                 % (link.init_node, link.term_node,
                    link.capacity, volume, link.linktype,
                    link.free_flow_time,link.B,link.power,link.length,
                    link.speed_limit,
                    link.change_id,
                     start_xy[0],
                    start_xy[1], end_xy[0], end_xy[1], 
                    linewidth,
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
    
    sys.stderr.write("Usage: " +progname + 
                     " [-s scaleFactor] [-zcm] [-u mod_filename]"
                     " netfilename nodefilename\n")

    sys.stderr.write("     -s ScaleFactor multiply co-ordines by ScaleFactor. Used as\n"
                     "                         e.g. -s 1000  for the Berlin TNTP data.\n"
                     "     -u mod_filename parse road upgrades from mod_filename\n"
                     "        and highlighed upgarded/added roads in red.\n"
                     "     -m do not read flow on stdin, just draw map\n"
                     "     -c  plot centroid connectors as aqua lines\n"
                     "     -z  plot location of zone centroids as purple dots\n"
                     "     -d  colour linktype 999 magenta for debugging reducenetwork.py\n")
    sys.exit(1)


def main():
    """
    main for nodenetflow2svg.py

    Usage: nodenetflow2svg.py [-s scaleFactor] [-u mods_file] [-zcmd] netfilename nodefilename 

      -s ScaleFactor multiply co-ordines by ScaleFactor. Used s
                     e.g. -s 1000  for the Berlin TNTP data.
     -m do not read flow on stdin, just draw map
     -c  plot centroid connectors as aqua lines
     -z  plot location of zone centroids as purple dots
     -u mods_file  - read upgarde modifications file mods_file and
                     draw map with upgarded/added roads highligted in red
     -d : colour linktype 999 magenta for debugging reducenetwork.py

      netfilename    is name of the net file defining node and links
      nodefilename   is name of node file giving x,y cooridinates for nodes

    Input on stdin is the flow output file giving volume on links on stdin,
    unless -m is specified then stdin is not read and map is drawn with no
    volume indicators.

    Output is SVG on stdout.

    Example usage:
    
    nodenetflow2svg.py SiouxFalls_net.txt SiouxFalls_node.txt < SiouxFalls_flow.txt
    
    """
    scaleFactor = None
    plot_zones = False
    plot_centroid_connectors = False
    read_flows = True
    highlight_upgrades = False
    debug999 = False

    try:
        opts,args = getopt.getopt(sys.argv[1:], "s:zcmu:d")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        if opt == "-s":  # scaelfactor
            scaleFactor = float(arg)
        elif opt == "-z": # plot zones
            plot_zones = True
        elif opt == "-c": # plot centroid connectors
            plot_centroid_connectors = True
        elif opt == "-m": # do not read flows, just draw map
            read_flows = False
        elif opt == "-u": # parse mods file ane plot upgardes in red
            highlight_upgrades = True
            mods_file = arg
        elif opt == "-d": # colour linktype 999 magenta for debugging
            debug999 = True
        else:
            usage(sys.argv[0])

    if len(args) != 2:
        usage(os.path.basename(sys.argv[0]))

    netfilename = args[0]
    nodefilename = args[1]

    net = parse_net_file(netfilename)
    nodexydict = parse_node_file(nodefilename)
    if scaleFactor != None:
        for (nodeid, (x,y)) in nodexydict.iteritems():
            nodexydict[nodeid] = (x*scaleFactor, y*scaleFactor)

    if read_flows:
        flows = parse_flow_file(sys.stdin)
    else:
        flows = None

    if highlight_upgrades:
        netmod = parse_mod_file(mods_file)
    else:
        netmod = None

    write_map_svg(sys.stdout, net, nodexydict, flows, plot_zones,
                  plot_centroid_connectors, netmod, debug999)


if __name__ == "__main__":
    main()
    
