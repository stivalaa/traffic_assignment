#!/bin/bash

# pbs launching script:
 	
# run on Tesla C1060 (nv4)
#PBS -l nodes=nv4
#PBS -q gpu
#PBS -N gpu_sssp_tesla
#PBS -l walltime=88:0:0
#PBS -l file=100MB


uname -a >&2

module load cuda/4.0
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

#/usr/bin/time ./sssp rmat.txt

/usr/bin/time ./sssp rmat.txt.chicagoregionalflows
