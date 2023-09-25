#!/usr/bin/bash
#PBS -N getorg__pt.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

module load GetOrganelle/1.7.4-pre2-foss-2020b 

get_organelle_config.py --add embplant_pt 

get_organelle_from_reads.py -1 READ1.fastq.gz -2 READ_2.trim.fastq.gz -o <path to output dir> -t 8 -F embplant_pt

