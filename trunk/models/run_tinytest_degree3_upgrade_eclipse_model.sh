#!/bin/sh
#
# $Id: run_tinytest_degree3_upgrade_eclipse_model.sh 637 2011-08-31 04:49:57Z astivala $
#
# Run the trivial upgrade subset_degree3 example data with MiniZinc model using
# Eclipse. The sed command is needed to edit the flatzinc output in order
# to swap the two lines, the first of which assigns 
# TripleMinusSum_VHT = TripleMinuSumVHT_list and the second assigns to the
# TripleMinusSumVHT_list: if we don't swap these so assignment to list is
# first then Eclipse gives an error messages.
#

mzn2fzn --output-to-stdout upgrade_subset_degree3_float.mzn tinytest_degree3.dzn   |   sed -n '/TripleMinusSumVHT =/{h;n;G;p;d;};p' |  eclipse -e "flatzinc:fzn_run(fzn_ic)"

