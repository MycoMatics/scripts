#!/usr/bin/bash
#PBS -N remap_canu.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=05:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be


ml BWA/0.7.17-intel-2018b

bwa mem -t 10 -x ont2d canu_WGS_DH_hesp.contigs.fasta ../allreads.trimmed.fastq.gz > nanopore_mapping/canu_mapping.sam