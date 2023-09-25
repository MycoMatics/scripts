#!/usr/bin/bash
#PBS -N name.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -e stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=05:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

module load xxxx

cd $PBS_O_WORKDIR
