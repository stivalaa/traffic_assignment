#!/bin/bash

# pbs launching script:
 	
#PBS -l nodes=4:ppn=4
#PBS -l pmem=1GB
#PBS -N mpi_tap_fw
#PBS -l walltime=96:0:0


cd $PBS_O_WORKDIR
set CONV_RSH = ssh

module load openmpi-gcc

time mpirun -report-bindings -npersocket 2 -np 8 ./tap_frankwolfe_mpi -n 2 -r 0.00001199 ~/traffic_assignment/trunk/testdata/ChicagoRegional/ChicagoRegional_net.txt ~/traffic_assignment/trunk/testdata/ChicagoRegional/ChicagoRegional_trips.txt ChicagoRegional_mods.txt ChicagoRegional_flows_

