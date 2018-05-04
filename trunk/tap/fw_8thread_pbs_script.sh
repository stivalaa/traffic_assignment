#!/bin/bash

# pbs launching script:
 	
#PBS -l nodes=1:ppn=8
#PBS -l pmem=2GB
#PBS -N tap_FrankWolfe
#PBS -l walltime=93:0:0


cd $PBS_O_WORKDIR
set CONV_RSH = ssh

uname -a >&2
date >&2

ITERATIONS=10000

echo  "threads = 0" >&2
/usr/bin/time ./tap_frankwolfe -i $ITERATIONS  ../testdata/ChicagoRegional/ChicagoRegional_net.txt  ../testdata/ChicagoRegional/ChicagoRegional_trips.txt  > /dev/null 

for threads in 1 2 3 4 5 6 7 8
do
  echo  "threads = $threads" >&2
  /usr/bin/time ./tap_frankwolfe_pthread -i $ITERATIONS -n $threads ../testdata/ChicagoRegional/ChicagoRegional_net.txt  ../testdata/ChicagoRegional/ChicagoRegional_trips.txt  > /dev/null 
done



times


