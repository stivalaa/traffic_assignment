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
# $Id: mktimestab.sh 845 2011-11-10 03:38:56Z astivala $ 
#

cat <<EOF
\begin{tabular}{llr}
  \hline
  Platform & Algorithm & Time (s) \\\\
  \hline
EOF
grep average | sed 's/reuse_same_arrays/reuse-memory/g' | sed 's/_/ /g' | while read line
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
  algorithm=`echo "${line}" | cut -f1 | cut -d' ' -f2- `
  ms=`echo "${line}" | awk '{print $(NF-1)}'`
  seconds=`echo "$ms / 1000" | bc -l`
  printf '  %s & %s & %.3f \\\\\n' $platform  "$algorithm"  $seconds
done | sort -t\& -k3,3n
cat <<EOF
  \hline
\end{tabular}
EOF

