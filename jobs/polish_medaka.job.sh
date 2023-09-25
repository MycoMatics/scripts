#!/usr/bin/bash
#PBS -N medaka.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

ml medaka/1.4.3-foss-2020b


cd $VSC_SCRATCH_VO_USER/minion/wgs_danny/qc_assemblies/

medaka_consensus -o $VSC_SCRATCH_VO_USER/minion/wgs_danny/qc_assemblies/medaka_polish_flyeass_out \
      -i $VSC_SCRATCH_VO_USER/minion/wgs_danny/trimmed/hespb_trimmed_q10_2000.fastq.gz \
      -d $VSC_SCRATCH_VO_USER/minion/wgs_danny/assembly/flye/hespb_flye.fasta \
      -m r103_sup_g507 \
      -t 96  
