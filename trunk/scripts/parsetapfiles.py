###############################################################################
#
# parsetapfiles.py - functions to parse traffic assignment problem files
#
# File:    parsetapfiles.py
# Author:  Alex Stivala
# Created: March 2011
#
# $Id: parsetapfiles.py 790 2011-10-06 02:49:50Z astivala $
#
###############################################################################

"""
Functions to parse the test data _node, _net, (traffic assignment input)
and flow (traffic assignment output) data in the format from

http://www.bgu.ac.il/~bargera/tntp/

"""

import sys

#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------

class Link:
    """
    Link is the class representing each link from the net file
    NB node numbers are numbered from 1

    BPR function: travel_time(Q) = T_0 * (1 + alpha*(Q/Q_max)^beta)
    """
    def __init__(self):
        init_node      = None  # from node number
        term_node      = None  # to node number
        capacity       = None  # Q_max in BPR function
        length         = None
        free_flow_time = None  # T_0 in BPR function
        B              = None  # alpha in BPR function
        power          = None  # beta in BPR function
        speed_limit    = None
        toll           = None
        linktype      = None

        # fields set by traffic assignment (not parsed)
        volume         = None  # volume on this link
        cost           = None  # cost on this link
    
class Net:
    """
    Net is the class for containing the data in the net file
    """
    def __init__(self):
        # Metadata
        num_zones       = None  # number of zones
        num_nodes       = None  # number of nodes
        first_thru_node = None  # first node that is not a zone (NB starts at 1)
        num_links       = None  # number of links
        # link data
        links           = None  # list of Link objects

class LinkFlow:
    """
    LinkFlow is the class for cotaining data from flow file
    """
    def __init__(self):
        init_node      = None  # from node number
        term_node      = None  # to node number
        volume         = None  # volume on link from TAP output
        cost           = None  # cost on link from TAP output

MODTYPE_ADD    = 1
MODTYPE_CHANGE = 2

class NetMod:
    """
    NetMod is the class for containg data from the mods file represneting
    network upgrades/changes.
    Note that the may  be multiple changes with the same change_id;
    they are parts of the same upgrade (each with its own subcost adding
    up to the total cost of the upgrde)
    """
    def __init__(self):
        change_id      = None  # change identifier
        modtype        = None  # MODTYPE_ADD, MODTYPE_CHANGE
        mod_link       = None  # Link object with added/changed link data
        project_cost   = None  # cost of this (part of) the change

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def parse_net_file(netfilename):
    """
    Parse the network file as described at http://www.bgu.ac.il/~bargera/tntp/

    E.g.
    <NUMBER OF ZONES> 24
    <NUMBER OF NODES> 24
    <FIRST THRU NODE> 1
    <NUMBER OF LINKS> 76
    <END OF METADATA>
    
    ~ 	Init node 	Term node 	Capacity 	Length 	Free Flow Time 	B	Power	Speed limit 	Toll 	Type	;
	1	2	25900.20064	6	6	0.15	4	0	0

    Parameters:
        netfilename - filename of net file to read 

    Return value:
        class Net object with data from file
    """
    maxnodenum = 0
    net = Net()
    net.links = []
    for line in open(netfilename):
        if line[0] == "~" or len(line.lstrip().rstrip()) == 0:
            pass
        else:
            line = line.rstrip()
            if line[0] == "<": #metadata
                tag = line[1:].split(">")[0]
                val = line[1:].split(">")[1].lstrip().rstrip()
                if tag == "NUMBER OF ZONES":
                    net.num_zones = int(val)
                elif tag == "NUMBER OF NODES":
                    net.num_nodes = int(val)
                elif tag == "FIRST THRU NODE":
                    net.first_thru_node = int(val)
                elif tag == "NUMBER OF LINKS":
                    net.num_links = int(val)
                elif tag == "END OF METADATA":
                    pass
                else:
                    sys.stderr.write("warning: unknown metadata tag: " + tag + "\n")
            else:
                sline = line.lstrip().split("\t")
                link = Link()
                link.init_node = int(sline[0])
                link.term_node = int(sline[1])
                link.capacity = float(sline[2])
                link.length = float(sline[3])
                link.free_flow_time = float(sline[4])
                link.B = float(sline[5])
                link.power = float(sline[6])
                link.speed_limit = float(sline[7])
                link.toll = float(sline[8])
                link.linktype = int(sline[9])
                net.links.append(link)
                maxnodenum = max(link.init_node, link.term_node, maxnodenum)
                link.volume = 0.0
                link.cost   = link.free_flow_time

    if net.num_links != len(net.links):
        sys.stderr.write("error: num_links %d but %d links\n" % (net.num_links,
                         len(net.links)))

    if maxnodenum != net.num_nodes:
        sys.stderr.write("error: max node num is %d but num_nodes is %d\n" %
                         (maxnodenum, net.num_nodes))

    return net
    

def parse_node_file(nodefilename):
    """
    Parse the node file which gives X and Y co-oridinates for each node

    Parameters:
       nodefilename - filename of the node file

    Return value:
       dict { nodeid : ( x, y) }
    """
    xydict = {}
    for line in open(nodefilename):
        if not line[0].isdigit():
            continue
        sline = line.split()
        nodeid = int(sline[0])
        x = float(sline[1])
        y = float(sline[2])
        if xydict.has_key(nodeid):
            sys.stderr.write("warning: duplicate node %d\n", nodeid)
        xydict[nodeid] = (x,y)
    return xydict

def parse_flow_file(flow_fh):
    """
    Parse the flow file giving volume and cost from output from traffic assign
    
    Paraemters:
       flow_fh - open (read) filehandle of the
                  flow file output by traffic assignment program

    Return value:
       list of LinkFlow objects giving volume and cost on each link

       <NUMBER OF NODES>       25
       <NUMBER OF LINKS>       76
       <END OF METADATA>

    E.g.:

       ~       Tail    Head    :       Volume  Cost    ;
               1       2       :       4494.5008499645437041   6.0008161234622576785   ;
               1       3       :       8119.1900669362912595   4.0086912217072878661  

    """
    flows = []
    for line in flow_fh:
        if line[0] == "~" or len(line.lstrip().rstrip()) == 0:
            pass
        else:
            line = line.rstrip().lstrip()
            if line[0] == "<": #metadata
                tag = line[1:].split(">")[0]
                val = line[1:].split(">")[1].lstrip().rstrip()
                if tag == "NUMBER OF NODES":
                    num_nodes = int(val)
                elif tag == "NUMBER OF LINKS":
                    num_links = int(val)
                elif tag == "END OF METADATA":
                    pass
                else:
                    sys.stderr.write("warning: unknown metadata tag: " + tag + "\n")
            else:
                sline = [s.lstrip().rstrip() for s in line.lstrip().split("\t")]
                flow = LinkFlow()
                flow.init_node = int(sline[0])
                flow.term_node = int(sline[1])
                if sline[2] != ":":
                    sys.stderr.write("warning: bad line '%s'\n" % line)
                flow.volume = float(sline[3])
                flow.cost = float(sline[4].rstrip(';'))
                flows.append(flow)

    if num_links != len(flows):
        sys.stderr.write("error: num_links %d but %d link flows\n" %
                         (num_links, len(flowslinks)))

    return flows


def parse_trips_file(tripsfilename, ignore_nonzero_od=False):
    """
    Parse the origin-destination demand file (trips file) which gives
    the demand fro travel between origin and destination (zone) pairs
    as described at http://www.bgu.ac.il/~bargera/tnt/p

    Paraemeters:
       tripsfilename- filename of the trips file to parse
       ignore_nonzero_od - if True, ignore nonzero demand on orig=dest

    Return value:
       dict of dicts { origin : { destination : demand } }

       
    Format of file  is e.g.:

    <NUMBER OF ZONES> 24
    <TOTAL OD FLOW> 360600.0
    <END OF METADATA>


    Origin  1
        1 :      0.0;     2 :    100.0;     3 :    100.0;     4 :    500.0;     5 :    200.0;
        6 :    300.0;     7 :    500.0;     8 :    800.0;     9 :    500.0;    10 :   1300.0;

    """
    demands = {}
    origin = None
    for line in open(tripsfilename):
        if line[0] == "~" or len(line.lstrip().rstrip()) == 0:
            pass
        else:
            line = line.rstrip().lstrip()
            if line[0] == "<": #metadata
                tag = line[1:].split(">")[0]
                val = line[1:].split(">")[1].lstrip().rstrip()
                if tag == "NUMBER OF ZONES":
                    num_zones = int(val)
                elif tag == "TOTAL OD FLOW":
                    total_od_flow = float(val)
                elif tag == "END OF METADATA":
                    pass
                else:
                    sys.stderr.write("warning: unknown metadata tag: " + tag + "\n")
            else:
                if line.split()[0] == "Origin":
                    origin = int(line.split()[1])
                    if origin > num_zones:
                        sys.stderr.write("warning: origin %d > num_zones %d\n"
                                         % (origin, num_zones))
                else:
                    sline = line.split(";")
                    sline.pop() # remove empty due to ; terminating line
                    for destflow in sline:
                        dfsplit = destflow.split(":")
                        dest = int(dfsplit[0])
                        if dest > num_zones:
                            sys.stderr.write("warning: dest %d > num_zonse %d\n"
                                             % (dest, num_zones))
                        demand = float(dfsplit[1])
                        if not demands.has_key(origin):
                            demands[origin] = {}
                        if (origin == dest and demand != 0.0 and
                            not ignore_nonzero_od):
                            sys.stderr.write("warning nonzero demand for orig = dest = %d: set to 0\n" % (origin))
                            demand = 0.0
                        demands[origin][dest] = demand
    return demands


def net_to_linkdict(net):
    """
    Convert Net object from parse_net_file to dict { (from,to) : Link }
    mapping a link indexed by (from,to) tuple to Link data
    
    Parameters:
      net - Net object as returned by parse_net_file()
    
    Return value:
      dict {(from,to):Link} as described above
    """
    return dict( [ ((link.init_node, link.term_node), link )
                   for link in net.links ] )
        

def net_to_graph(net):
    """
    Convert Net object from parse_net_file to graph represented
    (as per dijkstra.py) as dict of dicts where G[v][w] for any v,w
    is cost of edge from v to w. Here v and w are just integers (node numbers).

    Parameters:
       net - Net object as returned by parse_net_file()

    Return value:
       graph (dict of dicts) as described above
    """
    netgraph = dict((i, {}) for i in xrange(1,net.num_nodes+1))
    for link in net.links:
        netgraph[link.init_node][link.term_node] = link.cost
    return netgraph



def flow_to_linkdict(flows):
    """
    Convert LinkFlow object from parse_flow_file to dict { (from,to) : LinkFlow }
    mapping a link indexed by (from,to) tuple to LinkFlow data
    
    Parameters:
       flows - list of LinkFlow objects as returned by parse_flow_file()
    
    Return value:
       dict {(from,to):LinkFlow} as described above
    """
    return dict( [ ((linkflow.init_node, linkflow.term_node), linkflow)
                   for linkflow in flows ] )


def demands_to_demanddict(demands):
    """
    Convert dict of dicts from parse_trips_file() to dict keyed by
    (origin, dest) tuple

    Parameters:
        demands -    dict of dicts { origin : { destination : demand } }

    Return value:
        dict { (origin, destination) : demand }
    """
    ddict = {}
    for (orig, demand_dict) in demands.iteritems():
        for (dest, routeflow) in demand_dict.iteritems():
            assert(not ddict.has_key((orig, dest)))
            ddict[(orig, dest)] = routeflow
    return ddict


def parse_mod_file(modfilename):
    """
    Parse the network modification file

    This has the format (one change per line):
    
     ChangeId ChangeType Initnode  Termnode  Capacity 	Length 	FreeFlowTime  B Power	Speedlimit   Toll Type	
 
    The ChangeId must be suitable to be put in a filename as flows output
     files have it appended. 
    Note that multiple entries with a single ChangeId are allowed:
    this descrbies multiple modifications considered as part of a single
    indivisible network improvement. This is often used for example
    to add bidirectional links for a new road.
 
    The ChangeType can be 'ADD' or 'CHANGE'
 
    Note (at least for now) nodes have to already exist (ie cannot add nodes,
    only links between existing nodes). TODO allow new nodes. 
    TODO allow link deletion.

    Parameters:
        modfilename - filename of net mod file to read 

    Return value:
        list of NetMod objects with data from file
    """
    netmods = []
    for line in open(modfilename):
        if line[0] == "~" or len(line.lstrip().rstrip()) == 0:
            pass
        else:
            mod = NetMod()
            line = line.rstrip()
            sline = line.lstrip().split("\t")
            mod.change_id = sline[0]
            changetype_str = sline[1]
            if changetype_str == "ADD":
                mod.modtype = MODTYPE_ADD
            elif changetype_str == "CHANGE":
                mod.modtype = MODTYPE_CHANGE
            else:
                sys.stderr.write("unknown modification type %s\n"% changetype_str)
                sys.exit(1)
            mod.mod_link = Link()
            mod.mod_link.init_node = int(sline[2])
            mod.mod_link.term_node = int(sline[3])
            mod.mod_link.capacity = float(sline[4])
            mod.mod_link.length = float(sline[5])
            mod.mod_link.free_flow_time = float(sline[6])
            mod.mod_link.B = float(sline[7])
            mod.mod_link.power = float(sline[8])
            mod.mod_link.speed_limit = float(sline[9])
            mod.mod_link.toll = float(sline[10])
            mod.mod_link.linktype = int(sline[11])
            mod.project_cost = float(sline[12])
            netmods.append(mod)
            
    return netmods

    
def parse_individual_vht_from_tap_err(individual_vht_file):
    """
    Parse the VHT (total cost as vehcile hours travlled) from the stderr
    output of tap_frankwolfe_mpi for each individual change in mods file
    (input to tap_rankwolfe_mpi)

    Parameters:
       individual_vht_file -   file of the  the output
                               (stderr of tap_frankwolfe_mpi) of the tap
                               model with mods_file as input (indivual 
                               changes)
 
    Return value:
        dict { id : vht} of VHT for each modification
    
    """
    vht_dict = {} # dict { id : vht }
    for line in open(individual_vht_file):
        sline = line.split()
        if len(sline) > 3 and sline[3] == "VHT":
            vht_dict[sline[2]] = float(sline[5])
    return vht_dict


def parse_pairwise_vht_from_tap_err(pariwise_vht_file):
    """
    Parse the VHT (total cost as vehcile hours travlled) from the stderr
    output of tap_frankwolfe_mpi -p for each pair of changes in mods file
    (input to tap_rankwolfe_mpi)

    Parameters:
       pairwise_vht_file -   file of the  the output
                              (stderr of tap_frankwolfe_mpi -p) of the tap
                               model with mods_file as input  ,
                             giving pairwise changes
 
    Return value:
          dict { (id1, id2) : pairwise_vht }
               always id1 < id2 lexicographically
    
    """
    pair_vht_dict = {} # dict { (id1, id2) : pairwise_vht }
                       # always id1 < id2 lexicographically
    for line in open(pariwise_vht_file):
        sline = line.split()
        if len(sline) > 3 and sline[3] == "VHT":
            id1_and_id2_str = sline[2]
            pairwise_vht = float(sline[5])
            andsplit = id1_and_id2_str.split("_and_")
            id1 = andsplit[0]
            id2 = andsplit[1]
            pair_vht_dict[tuple(sorted([id1,id2]))] = pairwise_vht
    
    return pair_vht_dict


def parse_subset_vht_from_tap_err(subset_vht_file):
    """
    Parse the VHT (total cost as vehcile hours travlled) from the
    stderr output of tap_frankwolfe_mpi -s for each subset (size > 2)
    of changes in mods file (input to tap_rankwolfe_mpi)

    The id strings are spearated by a '+' character in this format e.g.

      03_03_0101+03_95_0001+03_96_0024

    Parameters:
       subset_vht_file -   file of the  the output
                              (stderr of tap_frankwolfe_mpi -s) of the tap
                               model with mods_file as input  ,
                             giving subset (size>2) changes
 
    Return value:
          dict { set  : vht }
             set is a set of upgrade id strings
    
    """
    set_vht_dict = {} # dict { set : pairwise_vht }
    for line in open(subset_vht_file):
        sline = line.split()
        if len(sline) > 3 and sline[3] == "VHT":
            idstrings = sline[2]
            vht = float(sline[5])
            idset = set(idstrings.split("+"))
            set_vht_dict[frozenset(idset)] = vht
    
    return set_vht_dict


