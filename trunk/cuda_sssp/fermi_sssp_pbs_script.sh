#!/bin/bash

# pbs launching script:
 	
# run on Fermi architecture (Tesla M2070)

#PBS -q gpu
#PBS -l nodes=1:ppn=8
#PBS -N gpu_sssp_fermi
#PBS -l walltime=23:0:0

uname -a >&2

module load cuda/4.2
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

${HOME}/cuda/sdk/C/bin/linux/release/deviceQuery >&2

/usr/bin/time ./sssp rmat.txt.chicagoregional
#cuda-memcheck ./sssp rmat.txt.chicagoregional
#/usr/bin/time ./sssp 512.graph.rmat.txt

