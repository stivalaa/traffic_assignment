#!/usr/bin/env python
##############################################################################
#
# convert_eclipse_output_to_id_list.py - convert output from Eclipse
#                                        running MiniZinc model
#                                        upgrade_subset_float.mzn etc.
#                                        to list of project ids
#
# File:    convert_eclipse_output_to_id_list.py
# Author:  Alex Stivala
# Created: August 2011
#
# $Id: convert_eclipse_output_to_id_list.py 711 2011-09-19 02:06:18Z astivala $
#
#
##############################################################################

"""
Parse the output from Eclipse solutoin of MiniZinc model and
convert to list of project ids from the model. Eclipse does not process
the MiniZinc output statement and just gives list of upgardes as a set
of integers, we need to parse the upgrade id string list from the model
and use these intergers as indices into that list to conver tback to
project id strings.

Usage:
     convert_eclipse_output_to_id_list.py [-c] data.dzn
  
     data.dzn - the MiniZinc .dzn data file to parse upgradeName list from

     -c : also include total cost of list of project ids and VHT reduction(NPV)

     Input is Eclipse output of solving model.mzn on stdin.
     Output is comma-delimited list of project ids on stdout.


Example usage:

mzn2fzn --output-to-stdout ../models/upgrade_subset_float.mzn  ../models/berlin_center_float.dzn  | eclipse -e "flatzinc:fzn_run(fzn_ic)" | ../scripts/convert_eclipse_output_to_id_list.py berlin_center_float.dzn
"""

import sys,getopt

compute_costs = False

try:
    opts,args = getopt.getopt(sys.argv[1:], "c")
except:
    sys.stderr.write("usage: %s [-c] data.dzn\n", sys.argv[0])
    sys.exit(1)
for opt,arg in opts:
    if opt == "-c":  
        compute_costs = True

if len(args) != 1:
    sys.stderr.write("usage: %s [-c] data.dzn\n", sys.argv[0])
    sys.exit(1)

dznfilename = args[0]

for line in open(dznfilename):
    if line.split()[0] == "upgradeName":
        # format of lists in MinZinc is same as Python so we can use eval
        # after stripping smicolon
        namelist = eval(line.split('=')[1].rstrip()[:-1])

    if compute_costs and line.split()[0] == "cost":
        costlist = [int(c) for c in eval(line.split('=')[1].rstrip()[:-1])]

for line in sys.stdin:
    # Eclipse outputs e.g. "upgradeSet = {3,5,6,7,8};" 
    # where the integers are (1-based) indices intto the upgradeName list
    sline = line.split('=')
    if len(sline) > 1 and sline[0].rstrip().lstrip() == 'upgradeSet':
        # similar to list, can eval() this into a set after sripping emicolon
        index_set = eval(sline[1].rstrip()[:-1])

    # Eclipse output for objective function looks like:
    # totalBenefitDollarsPV = 2780870.7012862656__2780870.7012862675;
    if compute_costs and len(sline) > 1 and (sline[0].rstrip().lstrip() == 'totalBenefitDollarsPV'
         or sline[0].rstrip().lstrip() == 'totalBenefitVHT'):
        objective_value_interval_str = sline[1].lstrip().rstrip()
        objective_value = int(round(float(objective_value_interval_str.split('__')[0])))


index_list = sorted(list(index_set))
names = [namelist[i-1] for i in index_list]
if compute_costs:
    costs = [costlist[i-1] for i in index_list]
    total_cost = sum(costs)

sys.stdout.write(','.join(names) )
if compute_costs:
    sys.stdout.write(' ' + str(total_cost) + ' ' + str(objective_value))
sys.stdout.write('\n')


