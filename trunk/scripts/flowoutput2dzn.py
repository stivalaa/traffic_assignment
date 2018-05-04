#!/usr/bin/env python
##############################################################################
#
# flowoutput2dzn - read mods file and VHT from tap_frankwolfe_mpi stderr 
#                  and generate .dzn file for  input to Zinc model
#                  to compute optimal upgrde subbset
#
# File:    flowoutput2dzn.py
# Author:  Alex Stivala
# Created: May 2011
#
# $Id: flowoutput2dzn.py 706 2011-09-15 06:40:01Z astivala $
#
#
##############################################################################

"""
Parse the stderr from tap_frankwolfe_mpi giving VHT for each individual
upgrade and also pairs of upgrades (to get difference between sum 
of VHT improvement in 2 upgrades and the VHT improvements resulting
from both upgrades in combination) and write .dzn format for Zinc model 
to solve the quadratic program to select optimal subset within
fixed budget

Usage:
   flowoutput2dzn.py [-i] [-t timeperiod] [-s id_list] 
                     [-p previously_built_list]
                     [-g potential_upgrades_list]
                      [-u subsets_tap_stderr_file]
                      mods_file nochange_vht_file 
                       invidual_tap_stderr_file 
                       pairwise_tap_stderr_file budget

   -i   : round all numerical values to integer for use with
          upgrade_subset_int.zinc

   -p previously_built_id_list : comma-delimited list of change ids
                       that have been previously built, used for
                        computing VHT pairwise adjustment but not
                         included in budget or VHT benefit
                         Sets preivoiuslyBuiltSet in dzn for
                         upgrade_subset_pv.dzn
                         Use -p0 to specifiy empty set

   -g potential_upgrades_list: command-delimited list of change ids
                       to be consdiered as potential upgrades,
                       sets potentialUpgradeSet in dzn for
                       upgrade_subset_pv.dzn
                       Use -g0 to specify all upgrades


   -s id_list   :     instead of all change ids in data
                      use only the ids in id_list which
                      is a comma-delimited list of change ids e.g. 
                      ber03,ber04,ber05,ber06,ber06a
                    
   -u subsets_tap_stderr_file   subsets_tap_stderr_file  is the output
                               (stderr of tap_frankwolfe_mpi) of the tap
                               model with mods_file as input (subset (size>2)
                               changes). If specified, output not just
                               pairwise PairMinusSumVHT but also degree 3
                               interactions TripleMinusSumVHT

   -t timeperiod : add a t = timeperiod; to .dzn file for use with
                     upgrade_subset_pv.zinc model. also add VHT2DollarFactor
                     and interestRate

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

   budget                    is the total fixed budget to fit costs of 
                             selected upgrades within

Output is .dzn zinc data text on stdout.

Example usage:

   flowoutput2dzn.py -i ChicagoRegional_mods.txt ChicagoRegional_fw_out.txt
                       ChicagoRegional_mods_goliath.err
                       ChicagoRegional_pairwise.err 50000000

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

VHT2dollarFactor = 10.00 * 365 # dollars/VHT * days/year
interestRate = 0.040           # annual interest rate 
     
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
                     " [-s id_list] [-i] [-t]\n"
                     " [-p preiviously_built_list]\n"
                     " [-g potential_upgrades_list]\n"
                     " [-u subsets_tap_stderr_file]\n"
                     " mods_file nochange_vht_file invidual_tap_stderr_file "
                     "pairwise_tap_stderr_file budget\n")
    sys.stderr.write("   -i : round all values to integer\n")
    sys.stderr.write("   -s id_list : subsets of ids in id_list only\n")
    sys.stderr.write("   -t timeperiod: add time period, intersrate to dzn\n")
    sys.stderr.write("   -u substs_tap_stderr_file: data for triples also\n")
    sys.stderr.write("   -p preiouvlys_built_list: list of previously built upgrades\n")
    sys.stderr.write("   -g potential_upgrades_list: list of potenti upgrades to consider\n")
    sys.exit(1)

def main():
    """
    See usage message in module header block
    """
    use_integer = False
    timeperiod = None
    specified_subset = None
    subset_vht_file = None
    output_triple_data = False
    previously_built_set = None
    potential_upgrade_set = None
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "it:s:u:p:g:")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        if opt == "-i":  # round all numverical values to integer
            use_integer = True
        elif opt == "-s": # use onlythis set  of ids not all
            specified_subset = set(arg.split(','))
        elif opt == "-t": # add t = timeperiod to .dzn 
            timeperiod = int(arg)
        elif opt == "-u": # parse subsets results and output degree 3 adjustments
            subset_vht_file = arg
            output_triple_data = True
        elif opt == "-p": # set of preiouvlsy built change ids
            previously_built_set = set(arg.split(','))
        elif opt == '-g': # set of potential upgrades
            potential_upgrade_set = set(arg.split(','))
        else:
            usage(sys.argv[0])

    if len(args) != 5:
        usage(sys.argv[0])

    mods_file = args[0]
    nochange_vht_file = args[1]
    individual_vht_file = args[2]
    pariwise_vht_file = args[3]
    budget = (args[4])

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

    # sort everything by change id so lists/arrays are ordered properly
    # and line up right i.e. cost[i] is cost of upgradeName[i] etc.
    idlist = sorted(vht_dict.keys())

    if specified_subset:
        indexlist = [idlist.index(projid) for projid in sorted(list(specified_subset))]
        xlist = [idlist[i] for i in indexlist]
        idlist = [projid for projid in idlist if projid in specified_subset]
        assert xlist == idlist

    if previously_built_set:
        if len(previously_built_set) == 1 and list(previously_built_set)[0] == '0':
            empty_previously_built_set = True
        else:
            empty_previously_built_set = False
            previously_built_indexlist = [idlist.index(projid) for projid in sorted(list(previously_built_set))]
            previously_built_xlist = [idlist[i] for i in previously_built_indexlist]
            previously_built_idlist = [projid for projid in idlist if projid in previously_built_set]
            assert previously_built_xlist == previously_built_idlist
    else:
        previously_built_indexlist = []

    if potential_upgrade_set:
        if len(potential_upgrade_set) == 1 and list(potential_upgrade_set)[0] == '0':
            potential_upgrade_set = set(idlist)
        potential_upgrade_indexlist = [idlist.index(projid) for projid in sorted(list(potential_upgrade_set))]
        potential_upgrade_xlist = [idlist[i] for i in potential_upgrade_indexlist]
        potential_upgrade_idlist = [projid for projid in idlist if projid in potential_upgrade_set]
        assert potential_upgrade_xlist == potential_upgrade_idlist
        if previously_built_set:
            if len(previously_built_set.intersection(potential_upgrade_set)) > 0:
                sys.stderr.write('ERROR: intersection of potential upgrade and previously biult sets nonempty\n')
                sys.exit(1)
    else:
        potential_upgrade_indexlist = []

    # dict { (id1, id2) : pairwise_vht }, always id1 < id2 lexicographically
    pair_vht_dict = parse_pairwise_vht_from_tap_err(pariwise_vht_file)

    # convert  to dict { (id1, id2) : orig - pairwise_vht }
    pair_deltavht_dict = dict([( (id1,id2), orig_vht - pair_vht ) for
                          ( (id1,id2), pair_vht ) in pair_vht_dict.iteritems()])

    if output_triple_data:
        # dict { set : subset_vht }
        subset_vht_dict = parse_subset_vht_from_tap_err(subset_vht_file)
        # conver to dict { set : orig_vht - subset_vht }
        subset_deltavht_dict = dict([(subset, orig_vht - subset_vht) for
                                     (subset, subset_vht) in 
                                     subset_vht_dict.iteritems()])

        # add pairwise into the more general subset deltaVHT dictionary
        for (pair, deltavht) in pair_deltavht_dict.iteritems():
            subset_deltavht_dict[frozenset(pair)] = deltavht

        # each subset_minus_sum_vht[i] is dict { subset : subset_delta_vht - 
        #                                                  sum(deltaVHt in subset) }
        # NB including changes for all subsets of size i,  2 <= i <= n
        MAXSUBSETSIZE = 3  # only doing triples for now 
        subset_minus_sum_vht = {} # key by subset size
        for subsetsize in xrange(2, MAXSUBSETSIZE+1): 
            subset_minus_sum_vht[subsetsize] = {}  # key is the subset
            idset = set(idlist)
            assert len(idset) == len(idlist) 
            for subset_tuple in itertools.combinations(idset, subsetsize):
                subset = frozenset(subset_tuple)
                singlesum = sum([orig_vht - vht_dict[uid] for uid in subset])
                pairsum = singlesum
                smallersubsetsum = pairsum
                for i in xrange(2, subsetsize):
                    # confusing: 'subset' is the whole subset 
                    # 'subsubset' is the subset of size i of the 'subset'
                    for subsubset_tuple in itertools.combinations(subset, i):
                        subsubset = frozenset(subsubset_tuple)
                        assert len(subsubset) == i
                        smallersubsetsum += subset_minus_sum_vht[i][subsubset]

                assert not subset_minus_sum_vht[subsetsize].has_key(subset)
                subset_minus_sum_vht[subsetsize][subset] = (
                    subset_deltavht_dict[subset]  - smallersubsetsum )
        

    outfh = sys.stdout
    outfh.write('% Generated by: ' + ' '.join(sys.argv) + '\n')
    outfh.write('% On: ' + timestamp + '\n')

    if timeperiod != None:
        outfh.write('VHT2dollarFactor = %f ;\n' % VHT2dollarFactor)
        outfh.write('interestRate = %f ;\n' % interestRate)
        outfh.write('t = %d ;\n' % timeperiod);

    
    outfh.write('n = ' + str(len(idlist)) + ';\n')
    if use_integer:
        outfh.write('Budget = ' + str(int(round(int(budget)))) + ';\n')
    else:
        outfh.write('Budget = ' + str(float(budget)) + ';\n')

    # since projects may have multiple changes need to sum up costs for total
    costlist = [ sum(change_costs) for change_costs in [
            [mod.project_cost for mod in netmod if mod.change_id == projid]
            for projid in idlist ] ]
    if use_integer:
        costlist = map(lambda x:int(round(x)), costlist)

    # handy that python list str format has same syntax as zinc data,
    # except have to replace ' with "
    outfh.write('upgradeName = ' + re.sub("'",'"',str(idlist)) + ';\n')
    outfh.write('cost = ' + re.sub("'",'"',str(costlist))+';\n')
    vhtlist = ([orig_vht - vht_dict[changeid] 
                     for changeid in idlist])
    if use_integer:
        vhtlist = map(lambda x:int(round(x)), vhtlist)
    outfh.write('benefitVHT = ' + str(vhtlist) + ';\n')

    if previously_built_set:
        if empty_previously_built_set:
            outfh.write('previouslyBuiltSet = {};\n')
        else:
            previously_built_index_set = set(previously_built_indexlist)
            assert len(previously_built_index_set) == len(previously_built_indexlist)
            outfh.write('previouslyBuiltSet = {' + ','.join([str(x+1) for x in previously_built_index_set]) + '};\n')

    if potential_upgrade_set:
        potential_upgrade_index_set = set(potential_upgrade_indexlist)
        assert len(potential_upgrade_index_set) == len(potential_upgrade_indexlist)
        outfh.write('potentialUpgradeSet = {' + ','.join([str(x+1) for x in potential_upgrade_index_set]) + '};\n')

    outfh.write('PairMinusSumVHT = [')
    for i in idlist:
        outfh.write('| ')
        for j in idlist:
            if i == j:
                pairvht = vht_dict[i]
                pair_delta = orig_vht - pairvht
                sum_delta = orig_vht-vht_dict[i] + orig_vht-vht_dict[j]
                abs_error = pair_delta - sum_delta;
            else:
                try:
                    if  i < j:
                        pairvht = pair_vht_dict[(i,j)]
                    else:
                        pairvht = pair_vht_dict[(j,i)]
                    pair_delta = orig_vht - pairvht
                    sum_delta = orig_vht-vht_dict[i] + orig_vht-vht_dict[j]
                    abs_error = pair_delta - sum_delta;
                except KeyError:
#                    sys.stderr.write(
#                        'no pairwise for (%s, %s), set abs_error=0\n' % (i,j))
                    abs_error = 0
            pair_minus_sum_vht = abs_error

            if use_integer:
                outstr = str(int(round(pair_minus_sum_vht)))
            else:
                outstr = str(float(pair_minus_sum_vht))
            outfh.write(outstr + ', ')
#            print 'x',i,j,pair_delta, sum_delta,pair_minus_sum_vht
        outfh.write('\n')
    outfh.write('|];\n')
            
    if output_triple_data:
        outfh.write("\n")
        outfh.write(
            "TripleMinusSumVHT_list = %array3d(Upgrades,Upgrades,Upgrades);\n")
        outfh.write("[\n")
        for i in idlist:
            for j in idlist:
                for k in idlist:
                    if i == j or i == k or j == k:
                        adjvalue = 0
                    else:
                        subset = frozenset([i,j,k])
                        assert len(subset) == 3
                        adjvalue = subset_minus_sum_vht[3][subset]
                    if use_integer:
                        outstr = str(int(round(adjvalue)))
                    else:
                        outstr = str(float(adjvalue))
                    outfh.write(outstr + ', ')
                outfh.write("\n")
            outfh.write("\n")
        outfh.write("];\n")
        
   
if __name__ == "__main__":
    main()
