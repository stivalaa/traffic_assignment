#!/bin/sh
#
# Make CSV file with x,y co-ordiantes of each update's centroid from
# the data in the rtab file.
#
# $Id: make_xy_csv_file_from_rtabs.sh 589 2011-08-22 04:26:18Z astivala $
#

header=""

  rtabfile=compare_pairwise_linkflow_summary.rtab

  if [ -z "${header}" ]; then
    header="ChangeId,CentroidX,CentroidY"
    echo "${header}" 
  fi

  grep -v '^#' $rtabfile | awk 'FNR > 1 {printf "%s,%f,%f\n",$1,$11,$12; printf "%s,%f,%f\n",$2,$13,$14;}' | sort | uniq
 

