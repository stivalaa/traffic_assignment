#!/bin/sh
#
# File:    mktimestab.sh
# Author:  Alex Stivala
# Created: July 2011
#
# mktimestab.sh - make LaTeX table from average times in sssp output
#
# Input is stdout from sssp on stdin, output is to stdout
#
# $Id: mktimestab.sh 484 2011-07-14 04:16:36Z astivala $ 
#

cat <<EOF
\begin{tabular}{llr}
  \hline
  Platform & Algorithm & Time (s) \\\\
  \hline
EOF
grep average | sed 's/_/ /g' | sed 's/&//g' | while read line
do
  echo "${line}" | fgrep 'Dijkstra' >/dev/null 2>&1
  if [ $? -eq 0 ]; then
    dijkstra=1
  else
    dijkstra=0
  fi
  echo "${line}" | fgrep 'internal' >/dev/null 2>&1
  if [ $? -ne 0 -a $dijkstra -eq 1 ]; then
    continue
  fi
  platform=`echo "${line}" | awk '{print $1}'`
  algorithm=`echo "${line}" | awk '{print $2}'`
  if [ $algorithm == "CUDA" ]; then
    algorithm=`echo "{$line}" | awk '{print $3}'`
  fi
  echo "${line}" | fgrep ' LLL ' > /dev/null 2>&1
  if [ $? -eq 0 ]; then
    algorithm="${algorithm} LLL"
  fi
  ms=`echo "${line}" | awk '{print $(NF-1)}'`
  seconds=`echo "$ms / 1000" | bc -l`
  printf '  %s & %s & %.3f \\\\\n' $platform  "$algorithm"  $seconds
done | sort -t\& -k3,3n
cat <<EOF
  \hline
\end{tabular}
EOF

