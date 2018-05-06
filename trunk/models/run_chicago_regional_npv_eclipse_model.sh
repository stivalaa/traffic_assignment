#!/bin/sh

mzn2fzn --output-to-stdout ~/traffic_assignment/trunk/models/upgrade_subset_npv.mzn  ~/traffic_assignment/trunk/models/chicago_pairwise_1only_npv.dzn | eclipse -e "flatzinc:fzn_run(fzn_ic)"
