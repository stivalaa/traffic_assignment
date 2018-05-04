#!/bin/bash


uname -a >&2
date >&2


RUNS=10
ITERATIONS=1000

echo  "threads = 0" >&2
runNum=1
while [ $runNum -le $RUNS ];
do
    echo "runNum = $runNum" >&2
    /usr/bin/time ./tap_frankwolfe -i $ITERATIONS  ../testdata/ChicagoRegional/ChicagoRegional_net.txt  ../testdata/ChicagoRegional/ChicagoRegional_trips.txt  > /dev/null 
    runNum=`expr $runNum + 1`
done  

for threads in 1 2 3 4 5 6 7 8
do
    echo  "threads = $threads" >&2
    runNum=1
    while [ $runNum -le $RUNS ];
    do
       echo "runNum = $runNum" >&2
      /usr/bin/time ./tap_frankwolfe_pthread -i $ITERATIONS -n $threads ../testdata/ChicagoRegional/ChicagoRegional_net.txt  ../testdata/ChicagoRegional/ChicagoRegional_trips.txt  > /dev/null 
    runNum=`expr $runNum + 1`
    done
done 

times

