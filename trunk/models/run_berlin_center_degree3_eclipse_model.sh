#!/bin/sh

mzn2fzn  --output-to-stdout upgrade_subset_degree3_float.mzn berlin_center_degree3.dzn | sed -n '/TripleMinusSumVHT =/{h;n;G;p;d;};p' |  eclipse -e "flatzinc:fzn_run(fzn_ic)"

