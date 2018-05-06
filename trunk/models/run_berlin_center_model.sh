#!/bin/sh

./upgrade_subset_int berlin_center_int.dzn |  ../scripts/verify_upgrade_subset_result.py -i ../data/berlin_center/berlin_center_mods.txt ../data/berlin_center/berlin_center_fw.err ../data/berlin_center/berlin_center_individual.err ../data/berlin_center/berlin_center_eyeball_significant_pairwise.err

