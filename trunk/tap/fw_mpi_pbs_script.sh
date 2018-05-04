#!/bin/bash

# pbs launching script:
 	
#PBS -l nodes=3:ppn=4
#PBS -l pmem=1GB
#PBS -N mpi_tap_fw
#PBS -l walltime=1:0:0


cd $PBS_O_WORKDIR
set CONV_RSH = ssh

module load openmpi-gcc


time mpirun -np 3 -npersocket 1  ./tap_frankwolfe_mpi -n 4 ~/traffic_assignment/trunk/testdata/SiouxFalls/SiouxFalls_net.txt ~/traffic_assignment/trunk/testdata/SiouxFalls/SiouxFalls_trips.txt  SiouxFalls_mods.txt SiouxFalls_flows_



