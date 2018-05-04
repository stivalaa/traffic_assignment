#!/usr/bin/env python
##############################################################################
#
# convert_eclipse_output_schedule_to_id_schedule.py - convert output from Eclipse
#                                        running MiniZinc model
#                                        upgrade_schedule_float.mzn etc.
#                                        to schedule of project ids
#
# File:    convert_eclipse_output_schedule_to_id_schedule.py
# Author:  Alex Stivala
# Created: August 2011
#
# $Id: convert_eclipse_output_schedule_to_id_schedule.py 726 2011-09-21 01:37:47Z astivala $
#
#
##############################################################################

"""
Parse the output from Eclipse solutoin of MiniZinc model and
convert to schedule of project ids from the model. Eclipse does not process
the MiniZinc output statement and just gives list of upgardes as a set
of integers, we need to parse the upgrade id string list from the model
and use these intergers as indices into that list to conver tback to
project id strings.

The output format is similar to that used in greedy_heruristic.sh
so it can also be parsed by run_tap_schedule_verification.py and
schedule_output_to_tex.py

Usage:
     convert_eclipse_output_schedule_to_id_schedule.py  [-c] data.dzn

     -c : also include total cost of list of project ids and VHT reduction(NPV)
     data.dzn - the MiniZinc .dzn data file to parse upgradeName list from

     Input is Eclipse output of solving model.mzn on stdin.
     Output is schedule of upgrade ids on stdout.


Example usage:
   mzn2fzn --output-to-stdout upgrade_schedule_independent.mzn chicago_regional_indepedent_schedule.dzn | eclipse -e "flatzinc:fzn_run(fzn_ic)" | ../scripts/convert_eclipse_output_schedule_to_id_schedule.py -c chicago_regional_indepedent_schedule.dzn
"""

import sys,getopt

compute_costs = False

try:
    opts,args = getopt.getopt(sys.argv[1:], "c")
except:
    sys.stderr.write("usage: %s [-c] data.dzn\n"% sys.argv[0])
    sys.exit(1)
for opt,arg in opts:
    if opt == "-c":  
        compute_costs = True
    else:
        sys.stderr.write("usage: %s [-c] data.dzn\n"% sys.argv[0])
        sys.exit(1)


if len(args) != 1:
    sys.stderr.write("usage: %s [-c] data.dzn\n"% sys.argv[0])
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
    # Eclipse outputs e.g.:
    #totalBenefitDollarsNPV = 156077380.35824692__156077380.35824704;
    #upgradeSet = array1d(1..5, [{7},{3},{2,4},{5},{6,8}]);
    # where the integers are (1-based) indices intto the upgradeName list
    sline = line.split('=')
    if len(sline) > 1 and sline[0].rstrip().lstrip() == 'upgradeSet':
        # use eval to convert to list of sets, each index_set_list[t] is upgrade set for period t
        index_set_list = eval('['+line.split('[')[1].split(']')[0]+']')

    if compute_costs and len(sline) > 1 and (sline[0].rstrip().lstrip() == 'totalBenefitDollarsNPV'):
        objective_value_interval_str = sline[1].lstrip().rstrip()
        objective_value = int(round(float(objective_value_interval_str.split('__')[0])))


num_periods = len(index_set_list)
for period in xrange(num_periods):
    sys.stdout.write('Schedule at period %d: ' % (period  ))
    index_set = index_set_list[period]
    index_list = sorted(list(index_set))
    names = [namelist[i-1] for i in index_list]
    if compute_costs:
        costs = [costlist[i-1] for i in index_list]
        total_cost = sum(costs)

    sys.stdout.write(','.join(names) )
    if compute_costs:
        sys.stdout.write(' (' + str(total_cost/1000) + ') ')
    sys.stdout.write('\n')

if compute_costs:
    sys.stdout.write('Total schedule benefit PV: %d\n' % (objective_value/1000))
