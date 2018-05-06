#!/bin/bash

# pbs launching script:
 	
# run on GTX280 (nv1) [or nv2, which has 2 of them]
#PBS -l nodes=nv1
#PBS -q gpu
#PBS -N gpu_sssp_gtx280
#PBS -l walltime=88:0:0
#PBS -l file=100MB


uname -a >&2

module load cuda/4.0
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

#/usr/bin/time ./sssp rmat.txt
/usr/bin/time ./sssp rmat.txt.chicagoregionalflows

