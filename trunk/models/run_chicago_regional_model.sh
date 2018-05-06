#!/bin/sh

time ../models/upgrade_subset_int  chicago_pairwiwse_1only_int.dzn | ../scripts/verify_upgrade_subset_result.py -i ../tap/ChicagoRegional_mods.txt ../tap/results/ChicagoRegional_fw_out.txt ../tap/results/ChicagoRegional_mods_goliath.err ../tap/ChicagoRegional_only1pairwise.err 
