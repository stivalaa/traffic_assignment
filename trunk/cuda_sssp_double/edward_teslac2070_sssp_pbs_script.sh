#!/bin/bash

# pbs launching script:
 	
# run on Fermi architecture (Tesla C2070) on edward.hpc.unimelb.edu.au cluster

# reserve 1 node with 8 cores to get entire node to maek sure no clash
# with any other gpu users on the node

#PBS -q gpu
#PBS -l nodes=1:ppn=8
#PBS -N gpu_sssp_teslac2070
#PBS -l walltime=23:0:0

uname -a >&2

module load cuda
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

${HOME}/cuda/sdk/C/bin/linux/release/deviceQuery >&2

/usr/bin/time ./sssp rmat.txt.chicagoregionalflows
#cuda-memcheck ./sssp rmat.txt.chicagoregional
#/usr/bin/time ./sssp 512.graph.rmat.txt


