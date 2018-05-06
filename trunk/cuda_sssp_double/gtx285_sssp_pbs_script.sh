#!/bin/bash

# pbs launching script:
 	
# run on GTX285 (nv3)
#PBS -l nodes=nv3
#PBS -q gpu
#PBS -N gpu_sssp_gtx285
#PBS -l walltime=88:0:0
#PBS -l file=100MB


uname -a >&2

module load cuda/4.0
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

#/usr/bin/time ./sssp rmat.txt

/usr/bin/time ./sssp rmat.txt.chicagoregionalflows
