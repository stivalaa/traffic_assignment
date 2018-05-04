#!/usr/bin/env python
###############################################################################
#
# csv2taptrip.py.py - Convert demand data in CSV format to TAP trip file format
#
# File:    csv2taptrip.py.py
# Author:  Alex Stivala
# Created: March 2011
#
# $Id: csv2taptrip.py 431 2011-06-29 00:28:05Z astivala $
#
###############################################################################

"""
Convert demand data in CSV format to TAP trip file
file in the format from

 http://www.bgu.ac.il/~bargera/tntp/

for use in traffic assignment programs (TAP).

See usage in docstring for main()

"""


import sys,os
import csv

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# multiply demand figures by this
DEMAND_MULT = 1  # TODO - work out what these figures really mean

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def csv2taptrip(csv_fh, tripfile_fh, maxzones=None):
    """
    Convert CSV demand data to TAP trip file format

    Parameters:
      csvfilename - open (read) filehandle to read CSV demand data from
      tripfile_fh - open (write) filehandle to write trip data to
      maxzones - maximum number of zones to convert or None for all

    Return value:
       None.
       
    Note the first line of the CSV file is some sort of 'header'
    which just has the zone numbers 1,...,n
    and the first entry in each subsequent line is the corresponding
    zone number.
    """
    reader = csv.reader(csv_fh)
    header = reader.next()
    num_zones = int(header[-1])
    if num_zones != len(header) - 1:
        sys.stderr.write("last number in header is not length\n")
        exit (1)
    if maxzones:
        num_zones = maxzones

    total_flow = 0
    od_dict = {}
    linecount = 0
    for zonedata in reader:
        origin = int(zonedata[0])
        od_dict[origin] = {}
        colcount = 0
        for dest_str in header[1:]:
            dest = int(dest_str)
            demand = float(zonedata[dest])
            od_dict[origin][dest] = demand * DEMAND_MULT
            total_flow += demand * DEMAND_MULT
            colcount += 1
            if maxzones and colcount == maxzones:
                break
        linecount += 1
        if maxzones and linecount == maxzones:
            break

    tripfile_fh.write("<NUMBER OF ZONES> " + str(num_zones) + "\n")
    tripfile_fh.write("<TOTAL OD FLOW> " + str(total_flow) + "\n");
    tripfile_fh.write("<END OF METADATA>\n")
    tripfile_fh.write("\n")
    origins = sorted(od_dict.keys())
    for origin in origins:
        tripfile_fh.write("Origin " + str(origin) + "\n")
        dests = sorted(od_dict[origin].keys())
        count = 0
        for dest in dests:
            tripfile_fh.write("%ld : %f;" % (dest, od_dict[origin][dest]))
            count += 1
            if count < 5:
                tripfile_fh.write("     ")
            else:
                tripfile_fh.write("\n")
                count = 0
        tripfile_fh.write("\n\n")
        

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    Print usage message and exit
    """

    sys.stderr.write("Usage: " +progname + " [maxzones] < csvfile > tripfile\n")
    sys.exit(1)


def main():
    """
    main for csv2taptrip.py

    Usage: csv2taptrip.py [maxzones] < csvfile > tripfile
 
    Input is CSV file on stdin and output to stdout is TAP trip file format
    If maxzonse is specified, it is the maximum number of zones to
    convert (zones from 1..mazxones converted, subsequent discarded).

    Example usage:
        bunzip2  -c ~/g12apps/toll/trunk/DoT-Disk-310310/PVV_CMB_EXTERNAL_Base_2008AM.csv.bz2 | ./csv2taptrip.py 2253  > melbourne_trips.txt
    """
    maxzones = None
    if len(sys.argv) > 2:
        usage(os.path.basename(sys.argv[0]))
    elif len(sys.argv) == 2:
        maxzones = int(sys.argv[1])


    csv2taptrip(sys.stdin, sys.stdout, maxzones)


if __name__ == "__main__":
    main()

