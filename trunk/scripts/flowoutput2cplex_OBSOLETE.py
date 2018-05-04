#!/usr/bin/env python
##############################################################################
#
# flowoutput2cplex - read mods file and VHT from tap_frankwolfe_mpi stderr 
#                    and generate input to CPLEX to compute optimal 
#                    upgrde subbset
#
# File:    flowoutput2cplex.py
# Author:  Alex Stivala
# Created: May 2011
#
# $Id: flowoutput2cplex_OBSOLETE.py 336 2011-06-01 07:19:09Z astivala $
#
#
##############################################################################

"""
Parse the stderr from tap_frankwolfe_mpi giving VHT for each individual
upgrade and also pairs of upgrades (to get difference between sum 
of VHT improvement in 2 upgrades and the VHT improvements resulting
from both upgrades in combination) and write CPLEX interactive input
to solve the quadratic program to select optimal subset within
fixed budget

Usage:
   flowoutput2cplex.py mods_file nochange_vht_file invidual_tap_stderr_file 
                       pairwise_tap_stderr_file budget

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

Output is cplex interactive optimizer text on stdout.

Example usage:

   flowoutput2cplex.py ChicagoRegional_mods.txt ChicagoRegional_fw_out.txt
                       ChicagoRegional_mods_goliath.err
                       ChicagoRegional_nods_pairwise.err 50000000

"""

import sys
from parsetapfiles import parse_mod_file,NetMod
     
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
                     "pairwise_tap_stderr_file budget\n")
    sys.exit(1)

def main():
    """
    See usage message in module header block
    """
    if len(sys.argv) != 6:
        usage(sys.argv[0])

    mods_file = sys.argv[1]
    nochange_vht_file = sys.argv[2]
    individual_vht_file = sys.argv[3]
    pariwise_vht_file = sys.argv[4]
    budget = float(sys.argv[5])

    netmod = parse_mod_file(mods_file)

    for line in open(nochange_vht_file):
        if line[:10] == "total cost":
            orig_vht = float(line.split()[3])
            break
        
    vht_dict = {} # dict { id : vht }
    for line in open(individual_vht_file):
        sline = line.split()
        if len(sline) > 3 and sline[3] == "VHT":
            vht_dict[sline[2]] = float(sline[5])

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
    
    sys.stdout.write('enter upgradesubset\n')
    sys.stdout.write('maximize ')
    idlist = sorted(vht_dict.keys())
    for i in xrange(len(idlist)):
        for j in xrange(i, len(idlist)):
            sys.stderr.write("cannot enter QIP in CPLEX this way, too bad\n")

    

    
if __name__ == "__main__":
    main()
