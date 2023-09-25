#!/usr/bin/bash
#PBS -N canu.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=48:00:00
#PBS -l mem=100gb
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

ml canu/2.1.1-GCCcore-10.2.0


cd $PBS_O_WORKDIR

canu -p canu_hespb -d canu genomeSize=17m -nanopore-raw $VSC_SCRATCH_VO_USER/minion/wgs_danny/trimmed/hespb_trimmed_q10_2000.fastq.gz useGrid=false
