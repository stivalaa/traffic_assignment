#!/bin/sh
mzn2fzn --output-to-stdout ../models/upgrade_subset_npv.mzn  ../models/berlin_center_npv.dzn  | eclipse -e "flatzinc:fzn_run(fzn_ic)"

