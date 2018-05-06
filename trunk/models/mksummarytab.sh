#!/bin/sh
###############################################################################
#
# mksummarytab.sh - make TeX table summarizing NPV of different algorithms
#
# File:    greedy_heuristic.sh
# Author:  Alex Stivala
# Created: Septmeber 2011
#
# Write to stdout a LaTeX table summarizing the net present value resulting
# from the verifcation (run_tap_schedule_verification.py) of the upgrade
# schedules resulting from different algorithms
#
# $Id: mksummarytab.sh 801 2011-10-12 05:56:52Z astivala $
# 
###############################################################################


# These three lists all have to line up i.e. ith element in list is for algorithm i

# the schedule from shcedule_output_to_tex.py for each algorithm
#SCHEDULE_TEX_FILES="chicago_heuristic_results.tex  chicago_independent_npv_schedule_results.tex	chicago_independent_schedule_results.tex chicago_heuristic_npv_results.tex"
SCHEDULE_TEX_FILES="chicago_independent_npv_schedule_results.tex	chicago_heuristic_npv_results.tex"

# the output of run_tap_schedule_verificatino.py for each algorithm
#VERIFICATION_OUTPUT_FILES="verification/chicago_greedy_schedule_verify.out verification_indpendent_npv/tap_npv_schedule_verify.o1360711 verification_independent/tap_schedule_verify.o1358875 verification_npv/chicago_greedy_npv_schedule_verify.out"
VERIFICATION_OUTPUT_FILES="verification_indpendent_npv/tap_npv_schedule_verify.o1360711 verification_npv/chicago_greedy_npv_schedule_verify.out"

# algorithm descriptions (use _ for space)
#DESCRIPTIONS="Greedy_Heuristic Linear_Model_NPV Linear_Model Greedy_Heuristic_NPV"
DESCRIPTIONS="Linear_Model Greedy_Heuristic"

echo '\begin{tabular}{lr}'
echo '\hline'
echo 'Algorithm & NPV (\$000s) \\'
echo '\hline'

num_algorithms=`echo "$DESCRIPTIONS" | awk '{print NF}'`
i=1
while [ $i -le $num_algorithms ];
do
    schedule_tex_file=`echo "$SCHEDULE_TEX_FILES" | awk "{print \\$$i}"`
    verification_output_file=`echo "$VERIFICATION_OUTPUT_FILES" | awk "{print \\$$i}"`
    description=`echo "$DESCRIPTIONS" | awk "{print \\$$i}" | tr _ ' '`
    cost=`grep Expenditure $schedule_tex_file | cut -d'&' -f 7 | tr -d '\\\'`
    benefit_pv=`grep '^Total schedule benefit' $verification_output_file | awk '{print $NF}'`
    npv=`echo "($benefit_pv - $cost * 1000)/1000" | bc -l`
    printf '%s & %.0f \\\\\n' "$description" $npv
    i=`expr $i + 1`
done | sort -t'&' -k2,2nr

echo '\hline'
echo '\end{tabular}'
