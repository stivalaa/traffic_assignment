#!/usr/bin/env python
##############################################################################
#
# reducenetwork.py - Reduce network by simplyfing paths of degree 2 nodes
#
# File:    reducenetwork.py
# Author:  Alex Stivala
# Created: Septmeber 2011
#
# $Id: reducenetwork.py 811 2011-10-25 04:52:37Z astivala $
#
#
##############################################################################

"""
Simplify the network by replacing paths along degree 2 nodes with the
same attributes with a single link. This is useful because the network file
is likely to be derived from GIS data (ESRI etc.) and the many links are
just there for geographically accuracy which we do not require.

See e.g. Klunder & Post (2006), however there they are ignoring attributes
and replacing a whole path of degree 2 nodes with a single link ignoring
the changes in attributes. We want to keep the attributes, specifically:

LINKTYPE
B
POWER    (this and B are dervied from CAPINDEX, so implied by CAPINDEX equal)
CAPACITY

so that only links where these values are all equal are reduced to a
single link. The single link resulting also has the sum of the following
values of the constitutent links

LENGTH
FREE_FLOW_TIME

so that calculations of costs (times) in TAP are correct.

Usage:
   reduncenetwork.py netfilename nodefilename out_netfilename out_nodefilename out_linkmap_filename

      netfilename    is name of the net file defining node and links
      nodefilename   is name of node file giving x,y cooridinates for nodes
      out_netfilename is name of output net file for reduced netwirk
      out_nodefilename is name of output node file for recued network
      out_linkmapefilename is name of output file containing mapping of new single links to multiple links

   Output is simplified network no out_netfilename and corresponding nodes
   on out_nodefilename (WARNING: theese files are overwritten if they exist).
   Also the mapping between new (added) links and the multiple (now removed) links they correspdond
   to is written to out_linkmapfilename in the format (one mapping perline):

   newlink oldink1;oldlink2,...

   where newlink and the oldlinks are fromnode,tonode tuples (i.e. node numbers delimited by comma)
   and the list of comma-delimited oldlink tuples is delmited by semcolons e.g.

   32867,39893  32867,8378;8378,9893;9893,10983;10983,39893

   Note that these node numbers are all in the OLD (original node numbering), unlike the output
   net and node files which have them renumbered.

   Becaus nodes are removed, but have to be sequentially numbered,
   node (and link) numbers no longer correspond to the originals.

    A mapping from original node numbers to renumbered (after deleting nodes
    no longer used) so as to remain sequential nuode numbers is written to stdout in form
   
    oldnum newnum

    one per line.
   

Example usage:

   reducenetwork.py  melbourne_net.txt  melbourne_node.txt melbourne_reduced_net.txt melbourne_reduced_node.txt > melbourne_reduce_number_map.txt
     
"""

import sys
import getopt
from time import strftime, localtime
from parsetapfiles import Link,Net,parse_net_file,parse_node_file,net_to_linkdict
from count_link_common_paths import  net_to_graph # this version handles Berlin data
import pape

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# epsilon for compaing floating point link attribute values
EPS = 1e-06

#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------

def has_same_attributes(link1, link2):
    """
    Return True if the two links have the same attributes for our purposes,
    ie it is OK to merge them together into one link

    Parameters:
       link1 - Link object
       link2 - Link object
    Return value:
       True iff link1 and link2 have compatible attributes
    """
    return (link1.linktype == link2.linktype and
            abs(link1.B - link2.B) < EPS and
            abs(link1.power - link2.power) < EPS and
            abs(link1.capacity - link2.capacity) < EPS)

def combine_link_attributes(linklist):
    """
    Return a new Link object combining links in linklist, which must have
    the same attributes according to has_same_attributes(), the new link
    has those and the sum of the free flow time and length of links in linkslist

    Uses the init_node in the first element of list as init_node and term_node
    in last element of list as term_node

    Parameters:
       linklist - list of Link objects
    Return value
       new Link object with summed free flow time and length
    """
    # print 'xxx_',[(link.init_node,link.term_node) for link in linklist]
    # print 'xxx0',[link.linktype for link in linklist]
    # print 'xxx1',[link.B for link in linklist]
    # print 'xxx2',[link.power for link in linklist]
    # print 'xxx3',[link.capacity for link in linklist]
    assert all([has_same_attributes(linklist[0], link) for link in linklist[1:]])
    newlink = Link()
    newlink.init_node = linklist[0].init_node
    newlink.term_node = linklist[-1].term_node
    newlink.capacity = linklist[0].capacity
#    newlink.linktype = linklist[0].linktype
    newlink.linktype = 999 # FIXME  special inktype for testing
    newlink.length         = sum(link.length for link in linklist)
    newlink.free_flow_time = sum(link.free_flow_time for link in linklist) # T_0 in BPR function 
    newlink.B              = linklist[0].B  # alpha in BPR function
    newlink.power          = linklist[0].power  # beta in BPR function
    newlink.speed_limit    = linklist[0].speed_limit
    newlink.toll           = linklist[0].toll
    newlink.skipover_count = 0
    return newlink

def iterative_dfs(graph, start, path=[]):
  ''' iterative depth first search from start
  '''
  q=[start]
  while q:
    v=q.pop(0)
    if v not in path:
      path=path+[v]
      q=list(graph[v])+q
  return path


def get_shortest_path(pred, orig, dest):
    """
    Get the shortest path to dest given the predecessor list pred
    from single-source shortest path algorithm

    Parameters:
       pred - list of predecessor nodes from a sssp algorithm with orig 
              as source
       source - source node
       dest  - destinatino node

    Return value:
      list of nodes on shortest path to dest from orig
    """
    path = []
    v = dest
    while v != orig:
        if v == -1:
            return None
        path.append(v)
        v = pred[v]
    path.append(orig)
    path.reverse()
    return path

def  linkcmp(link1, link2):
    """
    comparison function for Link objects used to sort list of Link objects
    by 'from' node ascending and within that by 'to' node ascending
    """
    if link1.init_node <  link2.init_node:
        return -1
    elif link1.init_node > link2.init_node:
        return 1
    else:
        if link1.term_node < link2.term_node:
            return -1
        elif link1.term_node > link2.term_node:
            return 1
        else:
            return 0
        

def node_is_shortcuttable(node_i, path, netgraph, rev_netgraph):
    """
    Given a index node_i in path and path (list of node numbers) containing node
    and the two graph (dict of dicts) structures
    where rev_netgraph is just graph with all links reversed, return
    True if the node has links only to and from its sucessor and
    predecessor in the aupplied path.

    NB this is more compicated that just 'degree 2' since  in thse
    networks there are directed edges and a link a a typeicla
    road conissts of forwrad and backward link.

    node_i - index of node in path
    path - list of nodes
    netgraph - dict G[v][w] is cost of link from v to w
    rev_netgraph - netgraph with all links reversed
    """
    assert node_i > 0 and node_i < len(path)-1
    node = path[node_i]
    outnodes = netgraph[node].keys()
    innodes = rev_netgraph[node].keys()
    if len(outnodes) > 2 or len(innodes) > 2:
        return False
    return set.union(set(innodes), set(outnodes)) == set([path[node_i-1],path[node_i+1]])

#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname +
                     " netfilename nodefilename out_netfilename out_nodefilename out_linkmapfilename\n")
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

    if len(args) != 5:
        usage(sys.argv[0])


    net_filename = args[0]
    nodefilename = args[1]
    out_net_filename = args[2]
    out_node_filename = args[3]
    out_linkmap_filename = args[4]

    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())

    net = parse_net_file(net_filename)
    nodexydict = parse_node_file(nodefilename)

    # build internal packed adjancey list for shortest path C module
    # and set link weights to link free flow time
    sys.stderr.write('building internal data...')
    netgraph = net_to_graph(net)
#!    sys.stderr.write( '[debug]: ' + str(len(net.links)) +'\n')

    for link in net.links:
        link.cost = link.free_flow_time

    # links must be sorted by from node ascending and then to node to match
    # the C extensino module for srhotest paths
    net.links.sort(cmp=linkcmp)

    pape.build_adjlist(netgraph, net.num_nodes+1, net.first_thru_node)
    pape.update_weights([link.cost for link in net.links])
    sys.stderr.write('done\n')

    # build reversed graph (all links reversed from original)
    rev_netgraph = dict((i,{}) for i in xrange(1, net.num_nodes+1))
    for (fromnode, nodedict) in netgraph.iteritems():
        for (tonode, cost) in nodedict.iteritems():
            rev_netgraph[tonode][fromnode] = cost

    linkdict = net_to_linkdict(net)

    for link in net.links:
        link.skipover_count = 0 # number of composed edges that skip over this one
        

    # We will go through all shortest paths between zones and add
    # special edges that span the maximal portion of
    # a path where all the nodes in that path have edges only
    # to and from their successor and precessor nodes in the path
    # NB this is more compicated that just 'degree 2' since  in thse
    # networks there are directed edges and a link a a typeicla
    # road conissts of forwrad and backward link.
    # Also we insist on all links that are colaeseced into a shortcut having
    # the same attributes
    shortcuts_dict = {} # dict {(startnode, endnode):[skipped link tuples]} for composed links to add
    for orig in xrange(1, net.num_zones+1):
        sys.stderr.write("processing origin %d\n" % orig)
        pred = pape.pape(orig)
        for dest in xrange(1, net.num_zones+1):
            if orig == dest:
                continue
            pathnodes = get_shortest_path(pred, orig, dest)
            if pathnodes == None:
                sys.stderr.write('warning: no path from %d to %d\n' % (orig, dest))
                continue
            if len(pathnodes) < 3:
                continue
            i = 0
            while i < len(pathnodes)-1:
                startnode = pathnodes[i]
                j = i+1
                nodes_along_shortcut = pathnodes[i:j+1] # includes start and end nodes
                # convert to list of edge tuples that the shortcut skips over
                skipped_link_tuples = zip(nodes_along_shortcut, nodes_along_shortcut[1:])
                skipped_links = [linkdict[link_tuple] for link_tuple in skipped_link_tuples]
                while (j < len(pathnodes)-1 and
                       node_is_shortcuttable(j, pathnodes, netgraph,rev_netgraph) and
#                       all([has_same_attributes(skipped_links[0],x) for x in skipped_links[1:]]) ):
                       has_same_attributes(skipped_links[0], skipped_links[-1])):
                    j += 1
                    nodes_along_shortcut = pathnodes[i:j+1] # includes start and end nodes
                    # convert to list of edge tuples that the shortcut skips over
                    skipped_link_tuples = zip(nodes_along_shortcut, nodes_along_shortcut[1:])
                    skipped_links = [linkdict[link_tuple] for link_tuple in skipped_link_tuples]
#                    print 'ppp',[l.capacity for l in skipped_links]
                endnode = pathnodes[j]
                if not has_same_attributes(skipped_links[0], skipped_links[-1]):
                    j -= 1
                    endnode = pathnodes[j]
                    skipped_links.pop()
                    skipped_link_tuples.pop()
                    nodes_along_shortcut.pop()
                    
                if (j - i >= 2):
                    if ( not shortcuts_dict.has_key((startnode, endnode)) or
                         len(skipped_link_tuples) > len(shortcuts_dict[(startnode,endnode)]) ):
                        for link in skipped_links:
                            link.skipover_count += 1
                        shortcuts_dict[(startnode,endnode)] = list(skipped_link_tuples)
#                        print 'aaa',pathnodes
#                        print 'bbb',startnode,endnode
#                        print 'ccc',[(l.init_node,l.term_node,l.skipover_count) for l in skipped_links]
#                        print 'zzz',[l.capacity for l in skipped_links]
                        assert (skipped_link_tuples[0][0] == startnode and
                                skipped_link_tuples[-1][1] == endnode)
                        assert all([has_same_attributes(skipped_links[0],x) for x in skipped_links[1:]])
                i = j+1
#    print 'yyy',shortcuts

    # TODO FIXME this whole thing has become a bit of a mess, should rewrite it to not be so
    # complicated  / inefficient
    
    # remove any subsumed shortcut edges
    remove_list =[]
    shortcut_tuples = list(shortcuts_dict.iterkeys())
    for i in xrange(len(shortcut_tuples)):
        st1 = shortcut_tuples[i]
        sl1 = shortcuts_dict[st1]
        for j in xrange(i+1, len(shortcut_tuples)):
            st2 = shortcut_tuples[j]
            sl2 = shortcuts_dict[st2]
            if set(sl2).issubset(set(sl1)):
                remove_list.append(st2)
            elif set(sl1).issubset(set(sl2)):
                remove_list.append(st1)
#    print 'zzzzzz', remove_list
    for st in remove_list:
        shortcuts_dict.pop(st)

    
    if len(shortcuts_dict) == 0:
        sys.stderr.write("WARNING: no composed edges could be constructed\n")
    else:
        sys.stderr.write('found %d composed edges, max len %d avg len %f\n' %
                     (len(shortcuts_dict), max([len(x) for x in shortcuts_dict.itervalues()]),
                      sum([len(x) for x in shortcuts_dict.itervalues()])/float(len(shortcuts_dict))))

    # write th elink mapping file
    linkmap_fh = open(out_linkmap_filename, "w")
    linkmap_fh.write('~ Generated by: ' + ' '.join(sys.argv) + '\n')
    linkmap_fh.write('~ On: ' + timestamp + '\n')
    linkmap_fh.write('~newLink\toldlink1;oldlink2;oldlink3;...\n')
    for (shortcut_tuple, oldlink_tuple_list) in shortcuts_dict.iteritems():
        linkmap_fh.write("%s\t%s\n" %
                         (str(shortcut_tuple),
                          ';'.join([str(linktuple ) for linktuple in oldlink_tuple_list])))
    linkmap_fh.close()
                                       
            
    # now make list of the shortcut links and into the list
    composed_links = []
    for shortcut_tuple in shortcuts_dict.iterkeys():
#        print 'qqq',shortcuts_dict[shortcut_tuple]
        composed_links.append(combine_link_attributes(
                            [linkdict[link_tuple] for link_tuple in shortcuts_dict[shortcut_tuple]]))

    # add the composed links, make sure to sort again
    net.links += composed_links
    net.links.sort(cmp=linkcmp)

    # make list of links, removing the ones that are skipped over by the 'shortcut' edges
    reduced_links = [link for link in net.links if link.skipover_count == 0]

#    sys.stderr.write('FIXME not removing links for testing\n')
#    reduced_links = [link for link in net.links] # XXX FIXME not removing links for debuggin.g..
    
    # create dict { node : True } for nodes that have a link out AND
    # one for links in, no longer counting links that have been skipped over,
    # but adding in the new (shortcut / composed) edges
    num_nodes = net.num_nodes
    node_linksout = {}
    for link in reduced_links:
        node_linksout[link.init_node] = True
    if len(node_linksout) < num_nodes:
        sys.stderr.write(
          "of %d nodes, %d now have no links out\n" %
          (num_nodes, num_nodes - len(node_linksout)))
    node_linksin = {}
    for link in reduced_links:
        node_linksin[link.term_node] = True
    if len(node_linksin) < num_nodes:
        sys.stderr.write(
          "of %d nodes, %d now have no links in\n" %
          (num_nodes, num_nodes - len(node_linksin)))
            

    removed_nodes = [node for node in xrange(1, num_nodes+1) if
                     not node_linksout.has_key(node) and not node_linksin.has_key(node)]

    assert all([node >= net.first_thru_node for node in removed_nodes])
    num_reduced_nodes = num_nodes - len ( removed_nodes )
    sys.stderr.write("removing %d nodes with no links in or out\n" % len(removed_nodes))

    # But node numbering has to be sequential so build a dict
    # { node_num : new_node_num } mapping the current node numbers,
    # for nodes with links only, to a new sequential numbering
    renumber_dict = dict([(i,i) for i in xrange(1,net.num_zones+1)]) # zones unchanged
    removed_nodes_dict = dict([(node, True) for node in removed_nodes])
    seqnum = 1
    # for node in xrange(first_thru_node, num_zones+num_nodes+1):
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
    outfh.write("<NUMBER OF ZONES> %d\n" % net.num_zones)
    outfh.write("<FIRST THRU NODE> %d\n" % net.first_thru_node)
    outfh.write("<NUMBER OF NODES> %d\n" % num_reduced_nodes)
    outfh.write("<NUMBER OF LINKS> %d\n" % len(reduced_links))
    outfh.write("<END OF METADATA>\n")
    outfh.write('~ from \tto\tcapacity\tlength\tftime\tB\tpower\tspeed\ttoll\ttype\n')

    for link in reduced_links:
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
    for node in xrange(1,  num_reduced_nodes+1):
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
