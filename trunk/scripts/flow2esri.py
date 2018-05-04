#!/usr/bin/env python
###############################################################################
#
# flow2esri.py - put flow (volumes) from TAP back into ESRI format
#
# File:    flow2esri.py
# Author:  Alex Stivala
# Created: October 2011
#
# $Id: flow2esri.py 818 2011-10-27 03:42:32Z astivala $
#
###############################################################################

"""

Take flow output from TAP and add into ESRI file containing the road network.


Requires that spatialite (v2.4) is installed in order to create access
spatialiate db from ESRI shapefiles and extract data from it and write
another ESRI file.

"""

import sys,os,subprocess

from parsetapfiles import parse_flow_file,LinkFlow

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    Print usage message and exit
    """

    sys.stderr.write("Usage: " +progname + " road_shapefile_in road_shapeefile_out \n")
    sys.stderr.write("    do not put .shp etc suffix on shapefile use basename\n")
    sys.stderr.write("    WARNING: road_shapefile_out is overwritten if it exists\n")
    sys.exit(1)


def main():
    """
    main for flow2esri.py

    Usage: flow2esri.py road_shapefile_in road_shapefile_out

    road_shapefile_in -    the ESRI shapefile containing roads (links) data
    road_shapefile_out -   the ESRI shapefile to write with roads data and now volumes 

    Flow file (TAP output) is read from stdin.
    
    WARNING: road_shapefile_out is overwritten if it exists
    
    Note shapefiles have the names

    X.shp
    X.dbf
    X.shx

    where X is the shapefile name.

    Example usage:

    renumber_nodes.py melbourne_number_map.txt < Melbourne_expanded_flows.txt | flow2esri.py ~/g12apps/toll/trunk/DoT-Disk-310310/HWY_Base2008AM_V090414_polyline  melbourne_flows_polyline 


    """
    if len(sys.argv) != 3:
        usage(sys.argv[0])
        sys.exit(1)
    in_shapefile = sys.argv[1]
    out_shapefile = sys.argv[2]

    sql = subprocess.Popen("spatialite", shell=True,
                           stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
    sql.stdin.write(".loadshp " + in_shapefile + " roads UTF-8\n")
    sql.stdin.write("alter table roads add column volume double precision;\n")
    sql.stdin.write("update roads set volume = 0.0;\n")
    sql.stdin.write("alter table roads add column volratio double precision;\n")
    sql.stdin.write("update roads set volratio = 0.0;\n")
    sql.stdin.write("alter table roads add column cost double precision;\n")
    sql.stdin.write("update roads set cost = 0.0;\n")

    sql.stdin.write("create index ab_index on roads(a,b)\n;") #otherwise update below very slow
    flows = parse_flow_file(sys.stdin)

    nzcount =  0
    nodedict ={}
    count = 0
    for lf in flows:
        if lf.volume != 0:
            sql.stdin.write("update roads set volume = %.15f, volratio = %.15f/capacity, cost = %.15f where A == %d and B == %d;\n" %
                            (lf.volume, lf.volume, lf.cost, lf.init_node, lf.term_node))
            nzcount += 1
            nodedict[(lf.init_node,lf.term_node)] = lf.volume
    sys.stderr.write("nzcount = %d unique nzcount = %d\n" % (nzcount, len(nodedict)))
    assert(nzcount == len(nodedict))
    sql.stdin.write("select count(*) from roads where volume != 0;\n")
    sql.stdin.write(".dumpshp roads Geometry " + out_shapefile + " UTF-8 linestring\n")
    sql.stdin.close()
    sys.stderr.write(str(sql.stdout.readlines()))
    sys.stderr.write('\n')
    sql.stdout.close()
    
if __name__ == "__main__":
    main()

