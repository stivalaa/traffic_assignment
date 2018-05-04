#!/usr/bin/env python
##############################################################################
#
# verify_upgrade_subset_result - compute total VHT benefit from 
#                                resulting subset from upgrade_subset_int
#                                Zinc model to verify result is correct.
#
# File:    verify_upgrade_subset_result.py
# Author:  Alex Stivala
# Created: May 2011
#
# $Id: verify_upgrade_subset_result.py 598 2011-08-23 05:57:42Z astivala $
#
#
##############################################################################

"""
Parse the stderr from tap_frankwolfe_mpi giving VHT for each individual
upgrade and also pairs of upgrades (to get difference between sum 
of VHT improvement in 2 upgrades and the VHT improvements resulting
from both upgrades in combination) and compute the total VHT benefit
from the subset of upgrades from the output of upgrade_subset_result_int
Zinc model to verify result is correct.

Usage:
   verify_upgrade_subset_result.py [-i] mods_file nochange_vht_file 
                       invidual_tap_stderr_file 
                       pairwise_tap_stderr_file budget <  uprade_subset.out


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

   -i  : round values to integers

Input is output of runnig the upgrade_subset_int model on stdin.

Example usage:

   upgrade_subset_int chicago_int.dzn | 
   verify_upgrade_subset_result.py -i ChicagoRegional_mods.txt ChicagoRegional_fw_out.txt
                       ChicagoRegional_mods_goliath.err
                       ChicagoRegional_pairwise.err

"""

import sys
import re
import getopt
from time import strftime, localtime

from parsetapfiles import parse_mod_file,NetMod,\
                          parse_individual_vht_from_tap_err,\
                          parse_pairwise_vht_from_tap_err
     
#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

     
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
                     " mods_file nochange_vht_file invidual_tap_stderr_file "
                     "pairwise_tap_stderr_file\n")
    sys.exit(1)

def main():
    """
    See usage message in module header block
    """
    use_integer = False
    try:
        opts,args = getopt.getopt(sys.argv[1:], "i")
    except:
        usage(sys.argv[0])
    for opt,arg in opts:
        if opt == "-i":  # round all numverical values to integer
            use_integer = True
        else:
            usage(sys.argv[0])

    if len(args) != 4:
        usage(sys.argv[0])

    mods_file = args[0]
    nochange_vht_file = args[1]
    individual_vht_file = args[2]
    pariwise_vht_file = args[3]


    netmod = parse_mod_file(mods_file)
    orig_vht = None
    for line in open(nochange_vht_file):
        if  line[:10] == "total cost" :
            orig_vht = float(line.split()[-1])
            break

    if orig_vht == None:
        for line in open(nochange_vht_file):
            if "total cost" in line and line.rstrip()[-4:] != " ms)":
                orig_vht = float(line.split()[-1])
        
    vht_dict = parse_individual_vht_from_tap_err(individual_vht_file)
   
    # dict { (id1, id2) : pairwise_vht }, always id1 < id2 lexicographically
    pair_vht_dict = parse_pairwise_vht_from_tap_err(pariwise_vht_file)

    # since projects may have multiple changes need to sum up costs for total
    costdict = {}
    for projid in vht_dict.iterkeys():
        costdict[projid] = sum([change_cost for change_cost in
            [mod.project_cost for mod in netmod if mod.change_id == projid]])
    if use_integer:
        costdict = dict(map(lambda x:(x[0], int(round(x[1]))),
                            costdict.iteritems()))

#    print  'q',costdict,  list(vht_dict.iterkeys())

    # Output of the model looks like:
    #   totalBenefitVHT = 45589
    #   ["03_02_9005", "03_03_0101", "03_95_0001", "07_06_0014", "07_96_0013", "07_97_0055"]
    # so we can treat the list of upgrade ids just as if it is a Python
    # list since it has the same syntax
    for line in sys.stdin:
        if line.split()[0] == "totalBenefitVHT":
            model_vht_str = line.split()[2]
            if use_integer:
                model_vht = int(model_vht_str)
            else:
                model_vht = float(model_vht_str)
        elif line[0] == '[':
            upgrade_list = eval(line)
        elif line.split()[0] == "totalCost":
            model_cost_str = line.split()[2]
            if (use_integer):
                model_cost = int(model_cost_str)
            else:
                model_cost = float(model_cost_str)

    sys.stdout.write("upgrade_list = " + str(upgrade_list) + "\n")

    individual_vht_benefits = [orig_vht - vht_dict[upgrade_id]
                               for upgrade_id in upgrade_list]
    if use_integer:
        individual_vht_benefits = map(lambda x:int(round(x)), 
                                      individual_vht_benefits)

#    print 'a',vht_dict
#    print 'x',model_vht, individual_vht_benefits,sum(individual_vht_benefits)

    upgrade_list.sort() # since pair_vht_dict keys (i,j) must have i < j
    pairwise_vht_adjustment = []
    for i in xrange(len(upgrade_list)):
        for j in xrange(i+1, len(upgrade_list)):
            try:
                pair_delta = orig_vht - pair_vht_dict[(upgrade_list[i],
                                                       upgrade_list[j])]
            except KeyError:
                pair_delta = None
            sum_delta = (orig_vht - vht_dict[upgrade_list[i]] +
                         orig_vht - vht_dict[upgrade_list[j]])
            if pair_delta: 
                abs_error =  pair_delta - sum_delta
            else:
                abs_error = 0
            pairwise_vht_adjustment.append(abs_error)
    if use_integer:
        pairwise_vht_adjustment = map(lambda x:int(round(x)),
                                      pairwise_vht_adjustment)
    pairwise_vht_delta = sum(pairwise_vht_adjustment)
#    print 'y', pairwise_vht_adjustment,pairwise_vht_delta
    recomputed_vht =  sum(individual_vht_benefits) + pairwise_vht_delta
    sys.stdout.write('model_vht = %d, recomputed_vht = %d\n'%(model_vht,
                                                              recomputed_vht))
    if model_vht - recomputed_vht == 0:
        sys.stdout.write("OK\n")
    else:
        sys.stdout.write("error, delta = %f\n"%(recomputed_vht-model_vht))

    recomputed_cost = sum(costdict[u] for u in upgrade_list)
    sys.stdout.write('model_cost = %d, recomputed_cost = %d\n'%
                     (model_cost, recomputed_cost))
    if model_cost - recomputed_cost == 0:
        sys.stdout.write("OK\n")
    else:
        sys.stdout.write('cost error, delta = %f\n'%
                         (recomputed_cost - model_cost))

                                                              


                                         
                                       
                                       
    
    

    
if __name__ == "__main__":
    main()
