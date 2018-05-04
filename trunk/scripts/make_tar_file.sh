#!/bin/sh
#
# make tar file of code etc. for distribution
#
# run from home directory (or other checkout of traffic_asignment/)
#
tar -c -z -v --exclude=tap/results --exclude-vcs -f tap.tar.gz traffic_assignment/trunk/tap traffic_assignment/trunk/scripts traffic_assignment/trunk/models traffic_assignment/trunk/python_extension traffic_assignment/trunk/Matlab
