#!/bin/bash

# pbs launching script:
 	
#PBS -l nodes=4:ppn=4
#PBS -l pmem=1GB
#PBS -N mpi_tap_fw
#PBS -l walltime=900:0:0


cd $PBS_O_WORKDIR
set CONV_RSH = ssh

module load openmpi-gcc


time mpirun -report-bindings -npersocket 2 -np 8 ./tap_frankwolfe_mpi  -s -n 2 -r 0.000011994618106 ../testdata/ChicagoRegional/ChicagoRegional_net.txt ../testdata/ChicagoRegional/ChicagoRegional_trips.txt ChicagoRegional_mods.txt ChicagoRegional_flows_ >  ChicagoRegional_mods_subsets.out 2> ChicagoRegional_mods_subsets.err 


