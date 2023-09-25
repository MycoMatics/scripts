#!/usr/bin/bash
#PBS -N flye_WGS_DH_hesp.$PBS_JOBID
#PBS -l nodes=1:ppn=64
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

ml Flye/2.8.3-GCC-10.2.0


cd $PBS_O_WORKDIR

flye --nano-raw allreads.fastq.gz -g 50M -o flye -t 8 