#!/usr/bin/env python
###############################################################################
#
# count_link_common_paths.py - count shortest paths through a pair of links
#
# File:    count_link_common_paths.py
# Author:  Alex Stivala
# Created: March 2011
#
# $Id: count_link_common_paths.py 784 2011-10-05 04:31:45Z astivala $
# 
###############################################################################

"""
Count the shortest paths (and total flow through those paths) that pass
through both of a supplied pair of links.

See usage decription in comment block for main.

This implementation uses a python extension written in C
(from ../python_extension/) to use
d'Esopo-Pape algorithm to make it faster

"""
import sys,os
from time import strftime,localtime,clock

import copy

import pape
from parsetapfiles import *


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def net_to_graph(net):
    """
    Convert Net object from parse_net_file to graph represented
    (as per dijkstra.py) as dict of dicts where G[v][w] for any v,w
    is cost of edge from v to w. Here v and w are just integers (node numbers).

    Parameters:
       net (in/OUT) - Net object as returned by parse_net_file()
                    duplicate entries in the links list are removed

    Return value:
       graph (dict of dicts) as described above
    """
    sys.stderr.write('[debug]: net_to_graph edges = ' + str(len(net.links))+'\n')
    netgraph = dict((i, {}) for i in xrange(1,net.num_nodes+1))
    delete_links = {} # dict { (initnode,termnode):seen } to delete after first
    for link in net.links:
        if (netgraph.has_key(link.init_node) and
            netgraph[link.init_node].has_key(link.term_node)):
            sys.stderr.write('WARNING: duplicate link %d -> %d\n' % 
                             (link.init_node, link.term_node))
            sys.stderr.write('         using first link only\n')
            delete_links[(link.init_node, link.term_node)] = False
        else:
            netgraph[link.init_node][link.term_node] = link.cost
    
    # now rebuild net.links without duplicate links (this happend on
    # Berlin data)
    net_links_copy = []
    for link in net.links:
#        net_links_copy.append(copy.deepcopy(link))
        # neither copy.copy() nor copy.deepcopy() actually seem to work at all
        # whether on lists or the objects in the list, have to do it manually.
        copylink = Link()
        copylink.init_node = link.init_node
        copylink.term_node = link.term_node
        copylink.capacity = link.capacity
        copylink.length = link.length
        copylink.free_flow_time = link.free_flow_time
        copylink.B = link.B
        copylink.power = link.power
        copylink.speed_limit = link.speed_limit
        copylink.toll = link.toll
        copylink.linktype = link.linktype
        copylink.cost   = link.cost
        net_links_copy.append(copylink)
        
    net.links = []
    for link in net_links_copy:
        if (not delete_links.has_key((link.init_node, link.term_node))
            or not delete_links[(link.init_node, link.term_node)]):
            net.links.append(link)
        delete_links[(link.init_node,link.term_node)] = True

    if net.num_links != len(net.links):
        sys.stderr.write('WARNING: had %d links, now %d after deleting duplicates\n' % (len(net_links_copy), len(net.links)))
    net.num_links = len(net.links)
    
    return netgraph


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
        # if v == -1:
        #     break
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
        

#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + 
                     " netfilename  demandfilename netmodfilename\n")
    sys.exit(1)


def main():
    """
    main for count_link_common_paths.py

    Usage: count_link_common_paths.py netfilename demandfilename netmodfilename

      netfilename    is name of the net file defining node and links
      demandfilename is name of Origin-Destination demand file
      netmodfilename is name of network modifications file, the modified
                     links here are used pairwise as the pairs to count
                     common paths through

      Input on stdin is flows file (output of tap_frankwolfe_mpi etc.)
      Output is common paths and flows on links for change pairs on stdout,
      with header line suitable for use in R read.table( ,header=T).
      
    Example usage:

    
    """
    if len(sys.argv) != 4:
        usage(os.path.basename(sys.argv[0]))

    netfilename = sys.argv[1]
    demandfilename = sys.argv[2]
    mods_file =sys.argv[3]

    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
    sys.stdout.write("# Generated by: " + " ".join(sys.argv) + "\n")
    sys.stdout.write("# On: " + timestamp + "\n")
    sys.stdout.write('ChangeId1\tChangeId2\tPaths\tFlow\n')
    

    sys.stderr.write('parsing net data...')
    net = parse_net_file(netfilename)
    linkdict = net_to_linkdict(net)
    sys.stderr.write('done.\n')

    sys.stderr.write('parsing trips data...')
    demands = parse_trips_file(demandfilename,ignore_nonzero_od=True)
    sys.stderr.write('done.\n')
    sys.stderr.write('parsing flow data...')
    flows = parse_flow_file(sys.stdin)
    sys.stderr.write('done.\n')

    # parse mods file and apply all the modifications to the network file
    netmod = parse_mod_file(mods_file)
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
                    net.links.append(mod.mod_link)
                    break
        else:
            raise ValueError('unknown modtype ' + str(mod.modtype))

    # links must be sorted by from node ascending and then to node to match
    # the C extensino module for srhotest paths
    net.links.sort(cmp=linkcmp)

    # build internal packed adjancey list for shortest path C module
    # and set link weights to costs from parsed flow file
    sys.stderr.write('building internal data...')
    flow_dict = flow_to_linkdict(flows)
    for link in net.links:
        try:
            link.cost = flow_dict[(link.init_node,link.term_node)].cost
        except KeyError:
            sys.stderr.write('no "cost" attribute in Link %d -> %d, '
                             'probably '
                             'due to flow input not including added link. '
                             'Setting cost to free flow time %f\n' % 
                             (link.init_node,link.term_node,
                              link.free_flow_time))
            link.cost = link.free_flow_time
    netgraph = net_to_graph(net)
    sys.stderr.write( '[debug]: ' + str(len(net.links)) +'\n')
    pape.build_adjlist(netgraph, net.num_nodes+1, net.first_thru_node)
    pape.update_weights([link.cost for link in net.links])
    sys.stderr.write('done\n')

    common_path_count = 0
    common_flow = 0


    # since projects (changeids) can have multiple link changes (or adds)
    # we get a list of all unique changeids, and for each a liset of
    # the changed/added links. Then we do all pairs of projects,
    # counting a common shortest path for that pair if it goes through both
    # [any of the links in changeid1] AND [any of th elinks in changeid2]
    changeid_list = sorted(list(set([mod.change_id for mod in netmod])))
    # dict { changeid : [list of (from,to)] }
    changeid_link_tuple_dict= dict([(changeid,[]) for changeid in changeid_list])
    for mod in netmod:
        changeid_link_tuple_dict[mod.change_id].append((mod.mod_link.init_node,
                                                        mod.mod_link.term_node))

    # changeid_list is sorted so these tuple indices always changeid1 < changeid2
    common_path_count_dict = {} # dict {(changeid1,changeid2):common_path_count }
    common_flow_dict = {} # dict {(changeid1,changeid2) : common_flow}
    for i in xrange(len(changeid_list)):
        changeid1 = changeid_list[i]
        for j in xrange(i+1, len(changeid_list)):
            changeid2 = changeid_list[j]
            common_path_count_dict[(changeid1,changeid2)] = 0
            common_flow_dict[(changeid1,changeid2)] = 0

    assert (len(list(common_path_count_dict.iterkeys())) == 
            len(changeid_list)*(len(changeid_list)-1)/2)

    sys.stderr.write("finding shortest paths...")
    t0 = clock()

    for (orig, demand_dict) in demands.iteritems():
        pred = pape.pape(orig)
        for (dest, routeflow) in demand_dict.iteritems():
            if orig == dest or routeflow == 0.0:
                continue
            pathnodes = get_shortest_path(pred, orig, dest)
            if len(pathnodes) < 2:
                continue
            path = zip(pathnodes, pathnodes[1:]) #convert to list of edge tuples

            # for all pairs of change ids, get common paths between any
            # of changeid1's links and changeid2's links
            # FIXME make this efficient rather than linear searches in path
            for i in xrange(len(changeid_list)):
                changeid1 = changeid_list[i]
                for j in xrange(i+1, len(changeid_list)):
                    changeid2 = changeid_list[j]
                    for link1_tuple in changeid_link_tuple_dict[changeid1]:
                        for link2_tuple in changeid_link_tuple_dict[changeid2]:
                            if link1_tuple in path and link2_tuple in path:
                                common_path_count_dict[(changeid1,changeid2)]+= 1
                                common_flow_dict[(changeid1,changeid2)]+=routeflow

    t1 = clock()
    sys.stderr.write("done (%f s).\n" % (t1 - t0))

    for (changeid1, changeid2) in common_path_count_dict.iterkeys():
        sys.stdout.write("%s\t%s\t%d\t%f\n" % 
                         (changeid1, changeid2,
                          common_path_count_dict[(changeid1,changeid2)],
                          common_flow_dict[(changeid1,changeid2)]))


if __name__ == "__main__":
    main()
