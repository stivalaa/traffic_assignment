#!/bin/bash

# pbs launching script:
 	
#PBS -l nodes=1:ppn=4
#PBS -l pmem=1GB
#PBS -N tap_FrankWolfe
#PBS -l walltime=23:0:0

uname -a >&2

cd $PBS_O_WORKDIR
set CONV_RSH = ssh

/usr/bin/time ./tap_frankwolfe_pthread -n4 ~/traffic_assignment/trunk/testdata/ChicagoRegional/ChicagoRegional_net.txt  ~/traffic_assignment/trunk/testdata/ChicagoRegional/ChicagoRegional_trips.txt  > ChicagoRegional_flows.txt


times


