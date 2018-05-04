#!/bin/bash

# pbs launching script:
 	
#PBS -l nodes=4:ppn=4
#PBS -l pmem=1GB
#PBS -N mpi_tap_fw_warm
#PBS -l walltime=23:0:0


cd $PBS_O_WORKDIR
set CONV_RSH = ssh

module load openmpi-gcc

time mpirun -npersocket 1 -np 4 ./tap_frankwolfe_mpi -w ChicagoRegional_flows.txt -n 4 -r 0.00001199 ~/traffic_assignment/trunk/testdata/ChicagoRegional/ChicagoRegional_net.txt ~/traffic_assignment/trunk/testdata/ChicagoRegional/ChicagoRegional_trips.txt ChicagoRegional_mods.txt ChicagoRegional_flows_warmstart_

