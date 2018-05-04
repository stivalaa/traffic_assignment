#!/usr/bin/env python
##############################################################################
#
# compare_pairwise_with_subsets - read mods file and VHT from 
#                  tap_frankwolfe_mpi stderr individual, pairwise, subsets
#                  and compare total VHT changes for individual+pairwise adj.
#                  with the actual VHT changes for each subset as gold stadndard
#
# File:    compare_pairwise_with_subsets.py
# Author:  Alex Stivala
# Created: June 2011
#
# $Id: compare_pairwise_with_subsets.py 776 2011-10-03 04:27:19Z astivala $
#
#
##############################################################################

"""
Parse the stderr from tap_frankwolfe_mpi giving VHT for each individual
upgrade and also pairs of upgrades (to get difference between sum 
of VHT improvement in 2 upgrades and the VHT improvements resulting
from both upgrades in combination) and also parse stderr for 
tap_frankwolfe mpi -s (all subsets) to get VHt improvements for each
subset. The latter is the gold stnadrd we compare the value computed
for VHT change form just individual and pairwise.
Also compare subsets of size 3,4,.. etc. unless -p is specfieid to only
do size 2 (pairwise) subsets.


Usage:
   compare_pairwise_with_subsets.py [-p] [-s id_list]
                    [  -i signficiant_pair_tap_stderr_file ]
                       mods_file nochange_vht_file 
                       invidual_tap_stderr_file 
                       pairwise_tap_stderr_file 
                       subsets_tap_stderr file


  -i signficiant_pair_tap_stderr_file read the VHT for the pairs in
                                      significant_tap_stderr_file, these
                                      are the ones computed pairwise for those
                                       pairs considered to interact.
                                       
  -p                  only do pairwise adjustments, not larger subsets

  -s id_list          instead of subsets of all change ids in data
                      for the SumAdjDeltaVHTx and RelDiffx columns,
                      do subsets of the set of ids in id_list which
                      is a comma-delimited list of change ids e.g. 
                      ber03,ber04,ber05,ber06,ber06a
                    

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

   subsets_tap_stderr_file    is the output
                               (stderr of tap_frankwolfe_mpi) of the tap
                               model with mods_file as input (subset (size>2)
                               changes)


Output is text suitable for R read.table(,header=T,stringsAsFactors=F)
on stdout.

Example usage:

   compare_pairwise_with_subsets.py  ChicagoRegional_mods.txt ChicagoRegional_fw_out.txt
                       ChicagoRegional_mods_goliath.err
                       ChicagoRegional_pairwise.err 
                       ChicagoRegional_mods_subsets.err

"""

import sys
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
# Constants
#
#-----------------------------------------------------------------------------

MAXSUBSETSIZE = 9   # max size of subsets to consider
     
#-----------------------------------------------------------------------------
#
# main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    print usage msg and exit
    """
    sys.stderr.write("usage: " + progname +  " [-p] [-s id_list] [-i signfiicant_pairs_tap_stderr_file] "+
                     " mods_file nochange_vht_file invidual_tap_stderr_file "
                     "pairwise_tap_stderr_file subsets_tap_stderr_file\n")
    sys.stderr.write("   -i significant_pairs_tap_stderr_file: read significant pair VHT values \n")
    sys.stderr.write("   -p : pairwise only, no subsets\n")
    sys.stderr.write("   -s id_list : subsets of ids in id_list only\n")
    sys.exit(1)

def main():
    """
    See usage message in module header block
    """
    pairwise_only = False
    specified_subset = None
    significant_pairs_tap_stderr_file = None

    try:
        opts,args = getopt.getopt(sys.argv[1:], "ps:i:")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        if opt == "-i":  # read signficiant interacting pairs VHT from file
            significant_pairs_tap_stderr_file = arg
        elif opt == "-p": # pairwise only
            pairwise_only = True
        elif opt == "-s": # use only subsets of this set for subset adjustments
            specified_subset = set(arg.split(','))
        else:
            usage(sys.argv[0])

    if len(args) != 5:
        usage(sys.argv[0])

    mods_file = args[0]
    nochange_vht_file = args[1]
    individual_vht_file = args[2]
    pariwise_vht_file = args[3]
    subset_vht_file = args[4]

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

    # dict { set : subset_vht }
    subset_vht_dict = parse_subset_vht_from_tap_err(subset_vht_file)
    # conver tto dict { set : orig_vht - subset_vht }
    subset_deltavht_dict = dict([(subset, orig_vht - subset_vht) for
                                 (subset, subset_vht) in 
                                 subset_vht_dict.iteritems()])

    # add pairwise into the more general subset deltaVHT dictionary
    for (pair, deltavht) in pair_deltavht_dict.iteritems():
        subset_deltavht_dict[frozenset(pair)] = deltavht

    if  significant_pairs_tap_stderr_file != None:
        # signfiicantly interacting pairs VHT
        significant_pairs_dict = parse_pairwise_vht_from_tap_err( significant_pairs_tap_stderr_file)
        # # convert  to dict { (id1, id2) : orig - pairwise_vht }
        # significant_pair_deltavht_dict = dict([( (id1,id2), orig_vht - pair_vht ) for
        #                            ( (id1,id2), pair_vht ) in significant_pairs_dict.iteritems()])
        
        
    outfh = sys.stdout
    outfh.write('# Generated by: ' + ' '.join(sys.argv) + '\n')
    outfh.write('# On: ' + timestamp + '\n')
    outfh.write('subsetDeltaVHT\tSumAdjPairwiseDeltaVHT\tRelDiff\tSumDeltaVHT\tSumRelDiff\tSubset\t')
    for i in xrange(2, MAXSUBSETSIZE+1):
        outfh.write('SumAdjDeltaVHT%d\tRelDiff%d' % (i, i))
        if i < MAXSUBSETSIZE:
            outfh.write('\t')
        else:
            outfh.write('\n')

    idlist = sorted(vht_dict.keys())
    
    # dict {(id1,id2): float(id1andid2_deltavht - (id1_deltavht + id2_deltavht))}
    #   always id1 < id2 lexicographically
    pair_minus_sum_vht = {}
    for k in xrange(len(idlist)):
        for l in xrange(k+1, len(idlist)):
            i = idlist[k]
            j = idlist[l]
            assert(i < j) # idlist is sorted so (k < l) => (i < j)
            try:
                if significant_pairs_tap_stderr_file != None:
                    pairvht = significant_pairs_dict[(i,j)]
                else:
                    pairvht = pair_vht_dict[(i,j)]
                pair_delta = orig_vht - pairvht
                sum_delta = orig_vht - vht_dict[i] + orig_vht - vht_dict[j]
                abs_error = pair_delta - sum_delta;
            except KeyError:
                sys.stderr.write('no pairwise for (%s, %s), set abs_error=0\n'
                                 % (i,j))
                abs_error = 0
            pair_minus_sum_vht[(i,j)] = abs_error

    if not pairwise_only:
        # each subset_minus_sum_vht[i] is dict { subset : subset_delta_vht - 
        #                                                  sum(deltaVHt in subset) }
        # NB including changes for all subsets of size i,  2 <= i <= n
        subset_minus_sum_vht = {} # key by subset size
        for subsetsize in xrange(2, MAXSUBSETSIZE+1): 
            subset_minus_sum_vht[subsetsize] = {}  # key is the subset
            if specified_subset != None:
                idset = specified_subset
            else:
                idset = idlist
            for subset_tuple in itertools.combinations(idset, subsetsize):
                subset = frozenset(subset_tuple)
                singlesum = sum([orig_vht - vht_dict[uid] for uid in subset])
                pairsum = singlesum
                # # include the pairwise adjustments
                # for pair in itertools.combinations(subset, 2):
                #     pairsum += pair_minus_sum_vht[tuple(sorted(pair))]
                # and all adjustments for subsets of size i where 2 < i < subsetsize
                smallersubsetsum = pairsum
                for i in xrange(2, subsetsize):
                    # confusing: 'subset' is the whole subset 
                    # 'subsubset' is the subset of size i of the 'subset'
                    for subsubset_tuple in itertools.combinations(subset, i):
                        subsubset = frozenset(subsubset_tuple)
                        assert len(subsubset) == i
                        smallersubsetsum += subset_minus_sum_vht[i][subsubset]

                assert not subset_minus_sum_vht[subsetsize].has_key(subset)
                subset_deltavht_value = subset_deltavht_dict[subset]
                subset_minus_sum_vht[subsetsize][subset] = (
                    subset_deltavht_value  - smallersubsetsum )
        

    # sys.stderr.write(str(pair_minus_sum_vht) + '\n')    #XXX
    # sys.stderr.write(str(subset_minus_sum_vht[3] ) +'\n')#XXX


    for (subset,subset_deltavht) in subset_deltavht_dict.iteritems():

        sumdeltavht = sum([deltavht_dict[uid] for uid in list(subset)])
        adj_sumdeltavht = sumdeltavht

        subset_idlist = sorted(list(subset))
        # now make adjustment for each pairwise change in the subset
        # according to the dictionary of pairwise VHT deltas
        for k in xrange(len(subset_idlist)):
            for l in xrange(k+1, len(subset_idlist)):
                i = subset_idlist[k]
                j = subset_idlist[l]
                assert(i < j) # subset_idlist is sorted so (k < l) => (i < j)
                adj_sumdeltavht += pair_minus_sum_vht[(i,j)]
#        print 'x',subset,subset_deltavht,sumdeltavht

        if not pairwise_only:
            # amke adjustmenst for each n-wise
            # change in the subset acording
            # to the dictionary of n-wise VHT deltas
            # NB including changes for all subsets of size i,  2 <= i <= n
            adj2plus_sum_deltavht = {}
#            for subsetsize in xrange(2, min(MAXSUBSETSIZE, len(subset))+1):
            for subsetsize in xrange(2,MAXSUBSETSIZE+1):
#                sys.stderr.write('yyy ' + str(subsetsize)+'\n')
                if subsetsize > 2:
                    adj2plus_sum_deltavht[subsetsize] = adj2plus_sum_deltavht[subsetsize-1]
                else:
                    adj2plus_sum_deltavht[subsetsize] = sumdeltavht
                for subset_tuple in itertools.combinations(subset, subsetsize):
                    nsubset = frozenset(subset_tuple)
                    assert len(nsubset) == subsetsize
                    try:
                        adj2plus_sum_deltavht[subsetsize] += subset_minus_sum_vht[subsetsize][nsubset]
                    except KeyError:
                        if specified_subset:
                            pass
                        else:
                            sys.stderr.write('missing data for subset ' + str(nsubset) + '\n')
                            raise KeyError

        outfh.write('%d\t%d\t%f\t%d\t%f\t%s\t' % (int(round(subset_deltavht)),
                                          int(round(adj_sumdeltavht)),
                                          (adj_sumdeltavht - subset_deltavht)/
                                          subset_deltavht,
                                          sumdeltavht,
                                          (sumdeltavht - subset_deltavht)/
                                          subset_deltavht,
                                          ','.join([str(x) for x in subset]) ))
        for subsetsize in xrange(2, MAXSUBSETSIZE+1):
            if not pairwise_only : #and subsetsize <= len(subset):
#                sys.stderr.write( 'xxx '+str(subset)+'\n')
                outfh.write('%d\t%f' % (
                        int(round((adj2plus_sum_deltavht[subsetsize]))),
                        (adj2plus_sum_deltavht[subsetsize] - subset_deltavht) /
                        subset_deltavht
                        ))
            else:
                outfh.write('NA\tNA')
            if subsetsize < MAXSUBSETSIZE:
                outfh.write('\t')
            else:
                outfh.write('\n')

    
if __name__ == "__main__":
    main()
