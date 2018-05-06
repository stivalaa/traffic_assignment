#!/bin/sh
mzn2fzn --output-to-stdout ../models/upgrade_subset_float.mzn  ../models/berlin_center_float.dzn  | eclipse -e "flatzinc:fzn_run(fzn_ic)"

