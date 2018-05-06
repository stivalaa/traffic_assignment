#!/bin/bash

# pbs launching script:
 	
# run on GTX 280 (nv1)
#PBS -l nodes=nv1
#PBS -q gpu
#PBS -N gpu_sssp_gtx280
#PBS -l walltime=1:0:0
#PBS -l file=100MB


uname -a >&2

module load cuda/2.3
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

/usr/bin/time ./sssp rmat.txt

