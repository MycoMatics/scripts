#!/usr/bin/bash
#PBS -N slaking_interactive.$PBS_JOBID
#PBS -l nodes=1:ppn=8
#PBS -l mem=26GB
#PBS -o slaking.stdout.$PBS_JOBID			
#PBS -e slaking.stderr.$PBS_JOBID
#PBS -l walltime=05:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be


source ~/.bashrc
conda activate pod5 &&
pod5 -h

