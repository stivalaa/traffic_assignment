#!/bin/sh

mzn2fzn --output-to-stdout ~/traffic_assignment/trunk/models/upgrade_subset_float.mzn  ~/traffic_assignment/trunk/models/chicago_pairwiwse_1only_float.dzn | eclipse -e "flatzinc:fzn_run(fzn_ic)"

