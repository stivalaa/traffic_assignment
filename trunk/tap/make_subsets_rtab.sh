#!/bin/sh
#
# Run the compare_pairwise_with_subset.py script and filter 
# to make table suitable
# for use with R read.table(,header=T) to analyse VHT computed from
# actual subsets and estimtated from only individual and pairwise.
#
# $Id: make_subsets_rtab.sh 557 2011-08-15 05:28:40Z astivala $
#
../scripts/compare_pairwise_with_subsets.py ChicagoRegional_mods.txt results/ChicagoRegional_fw_out.txt results/ChicagoRegional_mods_goliath.err results/ChicagoRegional_pairwise.err results/ChicagoRegional_mods_subsets.err | sed 's/froz/"froz/g' | sed 's/[)]/\)"/g'
