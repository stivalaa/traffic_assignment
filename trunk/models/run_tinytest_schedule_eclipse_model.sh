#!/bin/sh
#
# $Id: run_tinytest_schedule_eclipse_model.sh 616 2011-08-25 01:09:42Z astivala $
#
# Run the trivial upgrade schedule example data with MiniZinc model using
# Eclipse. The sed command is needed to edit the flatzinc output in order
# to swap the two lines, the first of which assigns 
# PairMinusSum_VHT = PairMinuSumVHT_list and the second assigns to the
# PairMinusSumVHT_list: if we don't swap these so assignment to list is
# first then Eclipse gives an error messages.
#

mzn2fzn --output-to-stdout upgrade_schedule_float.mzn tinytest_schedule.dzn   |   sed -n '/PairMinusSumVHT =/{h;n;G;p;d;};p' |  eclipse -e "flatzinc:fzn_run(fzn_ic)"

