#!/bin/bash

# pbs launching script:
 	
# run on Fermi architecture (Tesla C2050) (enrico)
# node enrico is 2.7GhZ quad core Nehalem server 12 GB RAM
# two Tesla C2050 GPUs

#PBS -q fermi
#PBS -N gpu_sssp_fermi
#PBS -l walltime=23:0:0
#PBS -l file=100MB

uname -a >&2

module load cuda/4.0
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

#${HOME}/cuda/sdk/bin/linux/release/deviceQuery >&2

/usr/bin/time ./sssp rmat.txt.chicagoregionalflows
#cuda-memcheck ./sssp rmat.txt.chicagoregional
#/usr/bin/time ./sssp 512.graph.rmat.txt

