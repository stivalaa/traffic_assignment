#!/bin/sh

if [ $# -ne 1 ]; then
  echo Usage: $0 pbs_job_number  >&2
  exit 1
fi

pbsjob=$1

while true;
do
# jobstatus does not use exit code for error, need to parse output
  jsout=`jobstatus m $pbsjob`
  echo `date` : $jsout
  echo "$jsout" | fgrep 'memory used by job'  >/dev/null 2>&1
  if [ $? -ne 0 ]; then
    break
  fi
  sleep 60
done
