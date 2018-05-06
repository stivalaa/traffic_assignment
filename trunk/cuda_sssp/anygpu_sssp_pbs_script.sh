#!/bin/bash

# pbs launching script:
 	
#PBS -l nodes=1
#PBS -q gpu
#PBS -N gpu_sssp_any
#PBS -l walltime=23:0:0
#PBS -l file=100MB


uname -a >&2

module load cuda/3.0
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

#cuda-memcheck ./sssp rmat.txt.chicagoregional
/usr/bin/time ./sssp rmat.txt.chicagoregional
#/usr/bin/time ./sssp ../cuda_apsp_withpaths/rmat.txt.winnipeg.unitcosts

