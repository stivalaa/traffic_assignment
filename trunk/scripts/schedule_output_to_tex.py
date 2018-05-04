#!/usr/bin/env python
##############################################################################
#
# schedule_output_to_tex.py - convert greedy_heurstic.sh output to LaTeX table
#
# File:    schedule_output_to_tex.py
# Author:  Alex Stivala
# Created: September 2011
#
# $Id: schedule_output_to_tex.py 739 2011-09-22 00:38:11Z astivala $
#
#
##############################################################################

"""
  Convert output from greedy_heuristic.sh like:

  ./greedy_heuristic.sh
  budgets:  1000 4000 1500 3000 5000
  total budget:  14500
  Set of changes to schedule:  03_03_0101,03_95_0001,07_06_0014,07_94_0027,07_96_0013,07_97_0055
  Total cost:  10384
  Total benefit PV:  158536
  Schedule at period 0: 07_96_0013 (748) 61411
  Schedule at period 1: 03_95_0001 (4000) 76078
  Schedule at period 2: 07_06_0014,07_94_0027 (1172) 23763
  Schedule at period 3: 03_03_0101 (465) 2780
  Schedule at period 4: 07_97_0055 (3999) 10184
  All scheduled

 to LaTex table. Input is on stdin, output to stdout

 Usage : ./greedy_heuristic.sh | schedule_output_to_tex.py netmods_file 
"""

import sys,os,re
from parsetapfiles import parse_mod_file,NetMod

if len(sys.argv) > 2 or len(sys.argv  ) < 2:
    sys.stderr.write("usage: schedule_output_to_tex.py netmods_file\n")
    sys.exit(1)

mods_file = sys.argv[1]
netmod = parse_mod_file(mods_file)
# since projects (changeids) can have multiple link changes (or adds)
# we get a list of all unique changeids
changeid_list = sorted(list(set([mod.change_id for mod in netmod])))

# since projects may have multiple changes need to sum up costs for total
costlist = [ sum(change_costs) for change_costs in [
        [mod.project_cost for mod in netmod if mod.change_id == projid]
        for projid in changeid_list ] ]
costlist = map(lambda x:int(round(x)/1000), costlist)
cost_dict = dict(zip(changeid_list, costlist))

period_cost_dict = {} # dict {period : totalcost}
project_dict = {} # dict {changeid : period}
for line in sys.stdin:
    sline = line.split()
    if len(sline) > 2 and sline[0] == 'budgets:':
        cline = line.split(':')
        budgets = cline[1].lstrip().rstrip().split()
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
    for changeid in upgrades.split(','):
        assert not project_dict.has_key(changeid)
        project_dict[changeid] = period
        period_total_cost += cost_dict[changeid]
#    print 'XXX',period_total_cost,cost
    assert abs(period_total_cost-  cost) < 2
    period_cost_dict[period] = period_total_cost

sys.stdout.write('\\begin{tabular}{rlrrrrrr}\n')
sys.stdout.write('\\hline\n')
sys.stdout.write('\\multicolumn{2}{c}{Time period $t$} & 1 & 2 & 3 & 4 & 5 & Total  \\\\\n')
sys.stdout.write('\\multicolumn{2}{c}{Budget (\$000s)}  & ')
for t in xrange(len(budgets)):
    sys.stdout.write(budgets[t])
    sys.stdout.write(' & ')
sys.stdout.write('%d' %  (sum(map(lambda x:int(x), budgets))))
sys.stdout.write('\\\\\n')
sys.stdout.write('$i$ & Project Id  \\\\\n')
sys.stdout.write('\\hline\n')

for i in xrange(len(changeid_list)):
    changeid = changeid_list[i]
    sys.stdout.write(str(i+1) + ' & ' + re.sub('_','-',changeid) + ' & ')
    for period in xrange(1+max(period_cost_dict.keys())):
        if project_dict.has_key(changeid) and project_dict[changeid] == period:
            sys.stdout.write('X')
        else:
            sys.stdout.write(' ')
        if period < max(period_cost_dict.keys()):
            sys.stdout.write('    &    ')
    sys.stdout.write(' \\\\\n')
sys.stdout.write('\\hline\n')

sys.stdout.write('\\multicolumn{2}{c}{Expenditure (\$000s) } & ')
for period in xrange(1+max(period_cost_dict.keys())):
    sys.stdout.write(str(period_cost_dict[period]))
    sys.stdout.write(' & ')
sys.stdout.write('%d' % (sum(period_cost_dict.values())))
sys.stdout.write(' \\\\\n')
sys.stdout.write('\\hline\n')
sys.stdout.write('\\end{tabular}\n')

