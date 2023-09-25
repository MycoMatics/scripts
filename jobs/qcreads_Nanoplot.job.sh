#!/usr/bin/bash


#PBS -l walltime=05:00:00							
#PBS -N rpb2_allreads_nanoplot									
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l nodes=1:ppn=all							
#PBS -m abe
								


module load NanoPlot/1.33.0-intel-2020b

cd /data/gent/vo/000/gvo00058/vsc43352/minion_runs/rpb2_20210625

NanoPlot -t 8 --fastq ./rpb2_allreads.fq.gz --loglength --plots dot --N50 --title rpb2_allreads -o ./QC_allreads