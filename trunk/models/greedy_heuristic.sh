#!/bin/sh
###############################################################################
#
# greedy_heuristic.sh - Greedy heurisic for road upgrade scheduling
#
# File:    greedy_heuristic.sh
# Author:  Alex Stivala
# Created: August 2011
#
#
# Solve the upgrade subset problem for the end of the planning periods
# and then with that subset only schedule the upgrades at each time
# period in the planning period.
# 
# This script uses the output of solving TAP with the upgrades and
# pairwise (only signfiicant) for each time period (different O-D demands)
# that must alkready by be available from runing tap_frankworlfe_mpi etc.
# It generates the .dzn data files for the models with python scripts
# and runs the MiniZinc models with this data with Eclipse
#
# The scripts/ directory with various python and shell scripts must be
# in PATH, as must mzn2fzn (from MiniZinc) and eclipse executables.
#
# $Id: greedy_heuristic.sh 714 2011-09-19 05:40:05Z astivala $
# 
###############################################################################


###############################################################################
#
# Parameters
#
###############################################################################

# number of planning periods (do not include original baseline)
NUM_PERIODS=4

# budget in each period (start at original baseline year)
# these are at present value (NB interest rate and VHT to dollar factor
# are currently set in the flowoutput2dzn.py script) 
BUDGETS="1000000.00 4000000.00 1500000.00 3000000.00 5000000.00"

# road network modificatnos file as used by tap_frankwolfe_mpi
MODS_FILE=$HOME/traffic_assignment/trunk/tap/ChicagoRegional_mods.txt

# original baseilne VHT result
BASELINE_ORIG_FW_OUT=$HOME/traffic_assignment/trunk/tap/results/ChicagoRegional_fw_out.txt

# baseline indivial modification TAP results
BASELINE_INDIVIDUAL_FW_OUT=$HOME/traffic_assignment/trunk/tap/results/ChicagoRegional_mods_goliath.err

# baseline pairwise (signficiant only) TAP results
BASELINE_PAIRWISE_FW_OUT=$HOME/traffic_assignment/trunk/tap/ChicagoRegional_only1pairwise.err

# Directory containing period T TAP results, where T = 1,2,...NUM_PERIODS 
# is appended to the directory name e.g. if this is /tap_data then
# we use /tap_data1 /tap_data2 etc.
PERIOD_TAP_RESULTS_DIR=$HOME/traffic_assignment/trunk/data/ChicagoRegional_tripmod

# filename in PERIOD_TAP_RESULTS_DIR of original (no upgrades) TAP results
PERIOD_TAP_ORIG_FW_OUT=ChicagoRegional_fw.out.txt

# filename in PERIOD_TAP_RESULTS_DIR of individual upgrade TAP results
PERIOD_TAP_INDIVIDUAL_FW_OUT=ChicagoRegional_mods_individual.err

# filename in PERIOD_TAP_RESULTS_DIR of pairwise (signfiacnt) upgrade TAP results
PERIOD_TAP_PAIRWISE_FW_OUT=ChicagoRegional_mods_pairwise.err


#directory containing the MiniZinc models
MODELS_DIR=$HOME/traffic_assignment/trunk/models

# where to store intermediate files and results
# WARNING will overwrite files used below, pid included but could still clobber
WORKDIR=/var/tmp

###############################################################################
#
# functions
#
###############################################################################

# return the set difference U \ V
# where U and V are sets repsented as comma-delimited lists
function set_difference() {
    U=$1
    V=$2
    ufile=`mktemp`
    vfile=`mktemp`
    echo $U | tr , '\n' | sort > $ufile
    echo $V | tr , '\n' | sort > $vfile
    udiffv=`comm -23 $ufile $vfile | tr '\n' , | sed 's/[,]$//'`
    rm $ufile $vfile
    echo "$udiffv"
    return 0
}

# return the set union of U and V
# where U and V are sets repsented as comma-delimited lists
function set_union() {
    U=$1
    V=$2
    if [ -z "$U" ]; then
      echo $V
    elif [ -z "$V" ]; then
      echo $U
    else
      echo ${U},${V}
    fi
    return 0
}


###############################################################################
#
# Main
#
###############################################################################

FINALPERIOD_DZN_FILE=$WORKDIR/endperiod$$.dzn


#
# first solve subset problem for end of planning periods with sum
# of budgets in each period
#
total_budget=0
budgets000=""
for i in $BUDGETS
do
    total_budget=`echo "$total_budget + $i" | bc`
    bdiv1000=`echo "$i / 1000" | bc`
    budgets000="$budgets000 $bdiv1000"
done
t=$NUM_PERIODS
flowoutput2dzn.py -t $t -g 0 -p 0 $MODS_FILE  ${PERIOD_TAP_RESULTS_DIR}$t/$PERIOD_TAP_ORIG_FW_OUT ${PERIOD_TAP_RESULTS_DIR}$t/$PERIOD_TAP_INDIVIDUAL_FW_OUT ${PERIOD_TAP_RESULTS_DIR}$t/$PERIOD_TAP_PAIRWISE_FW_OUT $total_budget  > $FINALPERIOD_DZN_FILE
upgrade_set_and_cost=`mzn2fzn --output-to-stdout $MODELS_DIR/upgrade_subset_pv.mzn $FINALPERIOD_DZN_FILE | eclipse -e "flatzinc:fzn_run(fzn_ic)" | convert_eclipse_output_to_id_list.py -c $FINALPERIOD_DZN_FILE`
upgrade_set=`echo "$upgrade_set_and_cost" | awk '{print $1}'`
upgrade_total_cost=`echo "$upgrade_set_and_cost" | awk '{print $2}'`
upgrade_benefit_pv=`echo "$upgrade_set_and_cost" | awk '{print $3}'`

#rm $FINALPERIOD_DZN_FILE

echo "budgets: " $budgets000
echo "total budget: " `echo "$total_budget / 1000" | bc`
echo "Set of changes to schedule: " $upgrade_set
echo "Total cost: " `echo "$upgrade_total_cost / 1000" |  bc`
echo "Total subset benefit PV (at end): " `echo "$upgrade_benefit_pv / 1000" | bc `

#
# now in each time period, solve the presnet value subset problem for
# that time period with the subset of upgrades not yet scheduled
#

t=0
previously_built_set=""
total_scheduled_benefit_pv=0
for budget in $BUDGETS
do 
    echo -n "Schedule at period ${t}: "
    PV_DZN_FILE=$WORKDIR/pv.t${t}.$$.dzn
    if [ $t -eq 0 ]; then
        flowoutput2dzn.py  $MODS_FILE  $BASELINE_ORIG_FW_OUT $BASELINE_INDIVIDUAL_FW_OUT $BASELINE_PAIRWISE_FW_OUT $budget  > $PV_DZN_FILE
        scheduled_set_and_cost=`mzn2fzn --output-to-stdout $MODELS_DIR/upgrade_subset_float.mzn $PV_DZN_FILE | eclipse -e "flatzinc:fzn_run(fzn_ic)" | convert_eclipse_output_to_id_list.py -c $PV_DZN_FILE`
        scheduled_set=`echo "$scheduled_set_and_cost" | awk '{print $1}'`
        scheduled_total_cost=`echo "$scheduled_set_and_cost" | awk '{print $2}'`
        scheduled_benefit_vht=`echo "$scheduled_set_and_cost" | awk '{print $3}'`
        scheduled_benefit_pv=`echo "$scheduled_benefit_vht * 3650" | bc`  # FIXME shoudl use VHT2DollarFactor from dzn file instead of 3650 hardcoded
    else
        flowoutput2dzn.py -t $t -g $upgrade_set -p $previously_built_set $MODS_FILE ${PERIOD_TAP_RESULTS_DIR}$t/$PERIOD_TAP_ORIG_FW_OUT ${PERIOD_TAP_RESULTS_DIR}$t/$PERIOD_TAP_INDIVIDUAL_FW_OUT ${PERIOD_TAP_RESULTS_DIR}$t/$PERIOD_TAP_PAIRWISE_FW_OUT $budget  > $PV_DZN_FILE
        scheduled_set_and_cost=`mzn2fzn --output-to-stdout $MODELS_DIR/upgrade_subset_pv.mzn $PV_DZN_FILE | eclipse -e "flatzinc:fzn_run(fzn_ic)" | convert_eclipse_output_to_id_list.py -c $PV_DZN_FILE`
        scheduled_set=`echo "$scheduled_set_and_cost" | awk '{print $1}'`
        scheduled_total_cost=`echo "$scheduled_set_and_cost" | awk '{print $2}'`
        scheduled_benefit_pv=`echo "$scheduled_set_and_cost" | awk '{print $3}'`
    fi
    scheduled_benefit_pv=`echo "$scheduled_benefit_pv / 1000" | bc`
    echo $scheduled_set "(`expr ${scheduled_total_cost} / 1000`)" $scheduled_benefit_pv
    upgrade_set=`set_difference $upgrade_set $scheduled_set` 
    previously_built_set=`set_union $previously_built_set $scheduled_set`
    total_scheduled_benefit_pv=`echo "$scheduled_benefit_pv + $total_scheduled_benefit_pv" | bc`
    if [ -z "$upgrade_set" ]; then
        echo "All scheduled"
#        rm $PV_DZN_FILE
        break
    fi
    t=`expr $t + 1`
#    rm $PV_DZN_FILE
done
if [ ! -z "$upgrade_set" ]; then
    echo "NOT SCHEDULED: " $upgrade_set
fi
echo "Total schedule benefit PV: " $total_scheduled_benefit_pv
