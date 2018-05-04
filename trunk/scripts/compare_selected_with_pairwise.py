#!/usr/bin/env python
##############################################################################
#
# compare_selected_with_pairwise - read mods file and VHT from 
#                  tap_frankwolfe_mpi stderr individual, selected
#                  (assumed interacting) pairwise, and all pairwise
#                  and compare total VHT changes for individual +
#                  selected_pairwise adj.
#                  with the  VHT changes adjusted for all pairwise (for sets
#                  to large to use actual subsets) 
#
# File:    compare_selected_with_pairwise.py
# Author:  Alex Stivala
# Created: June 2011
#
# $Id: compare_selected_with_pairwise.py 689 2011-09-13 06:07:50Z astivala $
#
#
##############################################################################

"""
Parse the stderr from tap_frankwolfe_mpi giving VHT for each individual
upgrade and also all pairs of upgrades (to get difference between sum 
of VHT improvement in 2 upgrades and the VHT improvements resulting
from both upgrades in combination). Compare the VHT adjusted for
only the selected pairwise interactions with all pairwise to compare
accuracy of using only the selected (predicted interacting) pairwise
interacations with using all pairwise adjustments.


Usage:
   compare_selected_with_pairwise.py 
                       [-m max_subset_size]
                       mods_file nochange_vht_file 
                       individual_tap_stderr_file 
                       pairwise_tap_stderr_file   < signficiant_pairs_list

   -m max_subset_size: only measure accuracy on subsets up to (including)
                       max_subset_size

   mods_file          is the file containg description of road upgrades
                      (also used as input in tap_frankwolfe_mpi program)

   nochange_vht_file  is the output (stderr of tap_frankwolfe_mpi)
                      of the TAP model with no modifications

   invidivudal_tap_stderr_file   is the output
                               (stderr of tap_frankwolfe_mpi) of the tap
                               model with mods_file as input (indivual 
                               changes)

   pairwise_tap_stderr_file    is the output
                               (stderr of tap_frankwolfe_mpi) of the tap
                               model with mods_file as input (pairwise 
                               changes)


Input on stdin is list of upgrade id pairs considered to be signficiantly
interacting pairwise, one on each line delimtied by whitespace. E.g:

ber03  ber04
ber05  ber06

Output is text suitable for R read.table(,header=T,stringsAsFactors=F)
on stdout.

Example usage:


"""

import sys
import re
import getopt
from time import strftime, localtime
import itertools

from parsetapfiles import parse_mod_file,NetMod,\
                          parse_individual_vht_from_tap_err,\
                          parse_pairwise_vht_from_tap_err
     


#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

DELTA_VHT_THRESHOLD = 1000 # ignore abs(DeltaVHT ) < this threshold

#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname +  " [-m maxsubsetsize] "
                     " mods_file nochange_vht_file invidual_tap_stderr_file "
                     "pairwise_tap_stderr_file < signficant_pair_list\n")

    sys.exit(1)

def main():
    """
    See usage message in module header block
    """
    max_subset_size = None

    try:
        opts,args = getopt.getopt(sys.argv[1:], "m:")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        if opt == "-m":
            max_subset_size = int(arg)
        else:
            usage(sys.argv[0])
    
    if len(args) != 4:
        usage(sys.argv[0])

    mods_file = args[0]
    nochange_vht_file = args[1]
    individual_vht_file = args[2]
    pariwise_vht_file = args[3]


    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())

    netmod = parse_mod_file(mods_file)

    orig_vht = None
    for line in open(nochange_vht_file):
        if line[:10] == "total cost":
            orig_vht = float(line.split()[3])
            break
    if orig_vht == None:
        for line in open(nochange_vht_file):
            if "total cost" in line and line.rstrip()[-4:] != " ms)":
                orig_vht = float(line.split()[-1])
        
    vht_dict = parse_individual_vht_from_tap_err(individual_vht_file)
    # conver to dict { id : orig_vht - vht}
    deltavht_dict = dict([(uid, orig_vht - vht) for (uid, vht) in
                          vht_dict.iteritems()])
   
    # dict { (id1, id2) : pairwise_vht }, always id1 < id2 lexicographically
    pair_vht_dict = parse_pairwise_vht_from_tap_err(pariwise_vht_file)
    # convert  to dict { (id1, id2) : orig - pairwise_vht }
    pair_deltavht_dict = dict([( (id1,id2), orig_vht - pair_vht ) for
                          ( (id1,id2), pair_vht ) in pair_vht_dict.iteritems()])

    # get list of signficiant pairs from stdin into
    # dict { (id1, id2) : True } with always id1 < id2 lexicographically
    significant_pair_dict = {}
    for line in sys.stdin:
        sline = line.split()
        if len(sline) != 2:
            sys.stderr.write("ERROR: bad pair line on stdin: '%s%'\n" % line)
            continue
        pair = tuple(sorted(sline))
        if significant_pair_dict.has_key(pair):
            sys.stderr.write("WARNING: duplicate pair %s\n" % str(pair))
        significant_pair_dict[pair] = True
        
    outfh = sys.stdout
    outfh.write('# Generated by: ' + ' '.join(sys.argv) + '\n')
    outfh.write('# On: ' + timestamp + '\n')
    outfh.write('numPairsUsed\tSumAdjPairwiseDeltaVHT\tsumAdjSignficiantPairwiseDeltaVHT\tRelDiff\tSumDeltaVHT\tSumRelDiff\tSubset\n')

    idlist = sorted(vht_dict.keys())

    # dict {(id1,id2): float(id1andid2_deltavht - (id1_deltavht + id2_deltavht))}
    #   always id1 < id2 lexicographically
    pair_minus_sum_vht = {}
    for k in xrange(len(idlist)):
        for l in xrange(k+1, len(idlist)):
            i = idlist[k]
            j = idlist[l]
            assert(i < j) # idlist is sorted so (k < l) => (i < j)
            pairvht = pair_vht_dict[(i,j)]
            pair_delta = orig_vht - pairvht
            sum_delta = orig_vht - vht_dict[i] + orig_vht - vht_dict[j]
            abs_error = pair_delta - sum_delta;
            pair_minus_sum_vht[(i,j)] = abs_error



    if max_subset_size == None:
        max_subset_size = len(vht_dict.keys()) # all subsets

    for subset_size in xrange(1, max_subset_size+1):
        for subset in itertools.combinations(vht_dict.keys(), subset_size):
            sumdeltavht = sum([deltavht_dict[uid] for uid in list(subset)])
            if abs(int(sumdeltavht)) < DELTA_VHT_THRESHOLD:
                sys.stderr.write("skipped %s  as %d < %d\n" % 
                              (str(subset), abs(int(sumdeltavht)),
                               DELTA_VHT_THRESHOLD))
                continue
            adj_sumdeltavht = sumdeltavht
            significant_adj_sumdeltavht = sumdeltavht
            subset_idlist = sorted(list(subset))
            # now make adjustment for each pairwise change in the subset
            # according to the dictionary of pairwise VHT deltas
            num_pairs_used = 0
            for k in xrange(len(subset)):
                for l in xrange(k+1, len(subset)):
                    i = subset_idlist[k]
                    j = subset_idlist[l]
                    assert(i < j) # subset_idlist is sorted so (k < l) => (i < j)
                    adj_sumdeltavht += pair_minus_sum_vht[(i,j)]
                    if significant_pair_dict.has_key((i,j)):
                        significant_adj_sumdeltavht += pair_minus_sum_vht[(i,j)]
                        num_pairs_used += 1
            try:
                RelDiff = ( (significant_adj_sumdeltavht - adj_sumdeltavht) /
                            adj_sumdeltavht )
            except ZeroDivisionError:
                RelDiff = float("NaN")

            try:
                SumRelDiff = ( (sumdeltavht - adj_sumdeltavht)/
                               adj_sumdeltavht )
            except ZeroDivisionError:
                SumRelDiff = float("NaN")

            outfh.write('%d\t%d\t%d\t%f\t%d\t%f\t%s\n' %
                        (num_pairs_used,
                         int(round(adj_sumdeltavht)),
                         int(round(significant_adj_sumdeltavht)),
                         RelDiff,
                         sumdeltavht,
                         SumRelDiff,
                        ','.join([str(x) for x in subset]) ))

    
if __name__ == "__main__":
    main()
