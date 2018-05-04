#!/usr/bin/env python
###############################################################################
#
# run_tap_schedule_verification.py - Run TAPs according to upgrade schedule
#
# File:    run_tap_schedule_verification.py
# Author:  Alex Stivala
# Created: September 2011
#
# $Id: run_tap_schedule_verification.py 738 2011-09-22 00:21:17Z astivala $
# 
###############################################################################
"""
 Read the upgrade schedule output from greedy_heuristic.sh and solve
 the TAP (using tap_frankwolfe_mpi) or read results from already
 exisintg TAP output for the upgrades at each time period in order to
 computer the actual VHT benefits according to the schedule.


 The tripmod files, network mod files, etc. must be the ones used
 by greedy_heruristic.sh (or otehr genrator of schedule) for this to make
 sense. At least for now, these data file locations are harcoded in 
 the Parmaters section belwo for Chicago Regional

 The schedule is read from stdin.

"""


import sys,os
import os.path
import re
import getopt
from time import strftime, localtime
import itertools

from parsetapfiles import parse_mod_file,NetMod,\
                          parse_individual_vht_from_tap_err,\
                          parse_pairwise_vht_from_tap_err,\
                          parse_subset_vht_from_tap_err     

#-----------------------------------------------------------------------------
#
# Constants / Parameters
#
#-----------------------------------------------------------------------------

VHT2dollarFactor = 10.00 * 365 # dollars/VHT * days/year
interestRate = 0.040           # annual interest rate 
     
HOME=os.getenv('HOME')
# road network modificatnos file as used by tap_frankwolfe_mpi
MODS_FILE=HOME+'/traffic_assignment/trunk/tap/ChicagoRegional_mods.txt'

NETWORK_FILE = HOME + "/traffic_assignment/trunk/testdata/ChicagoRegional/ChicagoRegional_net.txt"
TRIPS_FILE = HOME + "/traffic_assignment/trunk/testdata/ChicagoRegional/ChicagoRegional_trips.txt"

# original baseilne VHT result
BASELINE_ORIG_FW_OUT=HOME+'/traffic_assignment/trunk/tap/results/ChicagoRegional_fw_out.txt'

# baseline indivial modification TAP results
BASELINE_INDIVIDUAL_FW_OUT=HOME+'/traffic_assignment/trunk/tap/results/ChicagoRegional_mods_goliath.err'

#b baseline  pairwise modificatino TAP results
BASELINE_PAIRWISE_FW_OUT=HOME+'/traffic_assignment/trunk/tap/results/ChicagoRegional_pairwise.err'

#b baseline  subsets modificatino TAP results
BASELINE_SUBSETS_FW_OUT=HOME+'/traffic_assignment/trunk/tap/results/ChicagoRegional_mods_subsets.err'

# Directory containing period T TAP results, where T = 1,2,...NUM_PERIODS 
# is appended to the directory name e.g. if this is /tap_data then
# we use /tap_data1 /tap_data2 etc.
PERIOD_TAP_RESULTS_DIR=HOME+'/traffic_assignment/trunk/data/ChicagoRegional_tripmod'

# trip modification file without PERIOD_TAP_RESULTS_DIR
PERIOD_TRIPMODS_FILE="ChicagoRegional_tripmods.txt"

# filename in PERIOD_TAP_RESULTS_DIR of original (no upgrades) TAP results
PERIOD_TAP_ORIG_FW_OUT="ChicagoRegional_fw.out.txt"

# filename in PERIOD_TAP_RESULTS_DIR of individual upgrade TAP results
PERIOD_TAP_INDIVIDUAL_FW_OUT="ChicagoRegional_mods_individual.err"

# filename in PERIOD_TAP_RESULTS_DIR of pairwise upgrade TAP results
PERIOD_TAP_PAIRWISE_FW_OUT="ChicagoRegional_mods_pairwise.err"

# filename in PERIOD_TAP_RESULTS_DIR of all subsets of modifications TAP results
PERIOD_TAP_PAIRWISE_FW_OUT="ChicagoRegional_mods_subsets.err"


# location to store tap_frankwolfe_mpi output and flow files run from here
RESULTS_DIR="."

# basename for flows files from tap_frankworlfe_mpi
FLOWS_OUTPUT_BASENAME=RESULTS_DIR+"/ChicagoRegional_schedule_verifcation_flows_"

# basename for tap_frankwolfe_mpi stdout+stderr
TAP_OUTPUT_BASENAME=RESULTS_DIR+"/ChicagoRegional_schedule_verification_tap_output_"

# number of threads to use in tap_frankwolfe_pthread tap_frankwolfe_mpi
NUM_THREADS = 8

# relative gap target
RELATIVE_GAP_TARGET = 0.000011994618106

#-----------------------------------------------------------------------------
#
# functions
#
#-----------------------------------------------------------------------------

def parse_schedule(schedule_fh):
    """
    Parse the schedule output from greedy_heuristic.sh and return list of
    modifications to run TAP with at each time period

    Parameteres:
       shcedule_fh - open (read) filehandle of gready_heuristic.sh output
    Return value
       dict of lists, each key in dict is time period (0,1,2,...n)
        each inner list is list of change ids that are scheduled in that period
        ie dict { period : [ugparde list] }
    """
    upgrade_dict = {} # dict { period : [ upgarde list ] }
    for line in schedule_fh:
        sline = line.split()
        if len(sline) < 2 or sline[0] != 'Schedule':
            continue
        sline = line.split(':')
        period = int(sline[0].split()[3])
        sline = sline[1].split()
        upgrades = sline[0]
        cost = int(sline[1][1:-1])
        try:
            benefit = int(sline[2])
        except IndexError:
            benefit = None
        period_total_cost = 0
        assert not upgrade_dict.has_key(period)
        upgrade_dict[period] = upgrades.split(',')
    return upgrade_dict

        
def parse_orig_vht_from_nochange_file(nochange_vht_file):
    """
    Parse the original (no modifications) VHT from tap_frankwolfe_pthread file

    Parameters:
       nochange_vht_file - stder from tap_frankwolfe_pthread
    Return value
       VHT for the network prased from the file
    """
    orig_vht = None
    for line in open(nochange_vht_file):
        if line[:10] == "total cost":
            orig_vht = float(line.split()[3])
            break
    if orig_vht == None:
        for line in open(nochange_vht_file):
            if "total cost" in line and line.rstrip()[-4:] != " ms)":
                orig_vht = float(line.split()[-1])
    return orig_vht
            
            
#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname + "\n")
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

    if len(args) != 0:
        usage(sys.argv[0])

    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())

    netmod = parse_mod_file(MODS_FILE)
    upgrade_dict = parse_schedule(sys.stdin)

#    print 'xxx',upgrade_dict

    vht_change_list = [] # each vht_change_list[t] is VHT change in period t
    num_periods = max(upgrade_dict.keys()) + 1
    for period in xrange(num_periods):
        if period == 0:
            prev_changes = []
            # we always assume we already have this output 
            orig_vht = parse_orig_vht_from_nochange_file(BASELINE_ORIG_FW_OUT)
            if len(upgrade_dict[period]) == 1:
                vht_dict = parse_individual_vht_from_tap_err(BASELINE_INDIVIDUAL_FW_OUT)
                vht_change = orig_vht - vht_dict[upgrade_dict[period][0]]
            elif len(upgrade_dict[period]) == 2:
                # dict { (id1, id2) : pairwise_vht }, always id1 < id2 lexicographically
                pair_vht_dict = parse_pairwise_vht_from_tap_err(BASELINE_PAIRWISE_FW_OUT)
                vht_change = orig_vht - pair_vht_dict[tuple(sorted(upgrade_dict[period]))]
            elif len(upgrade_dict[period]) >= 3:
                # dict { set : subset_vht }
                subset_vht_dict = parse_subset_vht_from_tap_err(BASELINE_SUBSETS_FW_OUT)
                vht_change = orig_vht - subset_vht_dict[frozenset(upgrade_dict[period])]
            else:
                # no modifications
                vht_change = 0
            sys.stderr.write("period %d VHT change %f\n" % (period, vht_change))
            vht_change_list.append(vht_change)
        else:
            # at each period, solve TAP with all preivously scheduled upgrades
            # as the orig (baseline) of the network, then maeaure the VHT improvement
            # from the modigications scheduled for this period
            prev_changes.extend(upgrade_dict[period-1])
            if len(prev_changes) == 0:
                sys.stderr.write("TODO handle empty prev_changes list\n")
                sys.exit(99)
            tap_command_base = "tap_frankwolfe_mpi " + "-n " + str(NUM_THREADS) +\
                " -r " + str( RELATIVE_GAP_TARGET) +\
                " -t " + PERIOD_TAP_RESULTS_DIR + str(period)+"/" + PERIOD_TRIPMODS_FILE
            orig_output_name =  TAP_OUTPUT_BASENAME + str(period) + "_orig.err" 
            if not os.path.exists(orig_output_name):
                tap_command_orig = tap_command_base + " " +\
                    "-l " + ",".join(prev_changes) +\
                    " " + NETWORK_FILE +\
                    " " + TRIPS_FILE + " " + MODS_FILE +  " " +\
                    FLOWS_OUTPUT_BASENAME + " > " + orig_output_name + " 2>&1"
                sys.stderr.write("Running (orig): " + tap_command_orig + "\n")
                os.system(tap_command_orig)
            else:
                sys.stderr.write("using existing orig file %s\n" % orig_output_name)

            if len(prev_changes) == 1:
                vht_dict = parse_individual_vht_from_tap_err(orig_output_name)
                orig_vht = vht_dict[prev_changes[0]]
            elif len(prev_changes) >= 2:
                # dict { set : subset_vht }
                subset_vht_dict = parse_subset_vht_from_tap_err(orig_output_name)
                orig_vht = subset_vht_dict[frozenset(prev_changes)]
            else:
                assert False # FIXME
#            print 'aaaaa',orig_vht

   
            new_output_name = TAP_OUTPUT_BASENAME + str(period) + "_upgrade.err"
            if not os.path.exists(new_output_name):
                tap_command_new = tap_command_base + " " +\
                    "-l " + ",".join(prev_changes + upgrade_dict[period]) +\
                    " " + NETWORK_FILE +\
                    " " + TRIPS_FILE + " " + MODS_FILE +  " " +\
                    FLOWS_OUTPUT_BASENAME + " > " + new_output_name + " 2>&1"
                sys.stderr.write("Running (upgrade):" + tap_command_new + "\n")
                os.system(tap_command_new)
            else:
                sys.stderr.write("using existing upgrade file %s\n" % new_output_name)

            # dict { set : subset_vht }
            subset_vht_dict = parse_subset_vht_from_tap_err(new_output_name)
#            print 'xxxxx',frozenset(prev_changes + upgrade_dict[period])
#            print 'yyyyy',subset_vht_dict[frozenset(prev_changes + upgrade_dict[period])]
            vht_change = orig_vht - subset_vht_dict[frozenset(prev_changes + upgrade_dict[period])]

            sys.stderr.write("period %d VHT change %f\n" % (period, vht_change))
            vht_change_list.append(vht_change)

    # now we have all the VHT changes for each period, calculate NPV
    period_pv_list = [vht_change_list[period]*VHT2dollarFactor/(1+interestRate)**period for period in xrange(num_periods)]
    sys.stdout.write(str(vht_change_list)+'\n')
    sys.stdout.write(str(period_pv_list)+'\n')
    schedule_pv = sum(period_pv_list)
    sys.stdout.write('Total schedule benefit PV: ' + str(schedule_pv)+'\n')



if __name__ == "__main__":
    main()
