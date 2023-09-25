#!/usr/bin/bash
#PBS -N spades_HybAss_herp.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be




ml SPAdes/3.15.3-GCC-11.2.0


spades.py -o $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_hybridassembly/out_spadesassembly \
  -1 $VSC_SCRATCH_VO_USER/genomes/WMKJ_SRX7128323/trimmeddata/out_trimmomatic/SRR10432278_1_trimmed_paired.fastq \
  -2 $VSC_SCRATCH_VO_USER/genomes/WMKJ_SRX7128323/trimmeddata/out_trimmomatic/SRR10432278_2_trimmed_paired.fastq \
  --nanopore $VSC_SCRATCH_VO_USER/minion/wgs_danny/trimmed/hespb_trimmed_q10_2000.fastq.gz -t 96
  