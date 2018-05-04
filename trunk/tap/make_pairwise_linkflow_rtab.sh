#!/bin/sh
#
# Make the R table input file compare_pairwise_linkflow_summary.rtab (stdout)
# used (possibly after converting to csv with 
# scripts/pairwisesummaryrtab2csv.sh)  to train model to predict which 
# changes have pairwise interactions signfinactly differnet form sum of
# indivudual VHT chnages
#
# $Id: make_pairwise_linkflow_rtab.sh 460 2011-07-06 05:16:17Z astivala $
#
./summarize.sh  | ../scripts/compare_pairwise_linkflow_summary.py  ../testdata/ChicagoRegional/ChicagoRegional_net.txt ../testdata/ChicagoRegional/ChicagoRegional_node.txt ChicagoRegional_mods.txt results/ChicagoRegional_flows.txt results/ChicagoRegional_flows_

