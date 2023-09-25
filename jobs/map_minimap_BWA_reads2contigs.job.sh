#!/usr/bin/bash
#PBS -N remap_canu.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be
      
  # MINIMAP2-SAMTOOLS: index draft genome and align basecalled reads to draft
    ml minimap2/2.24-GCCcore-11.2.0
    ml SAMtools/1.14-GCC-11.2.0
    cd $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_hybridassembly
    minimap2 -ax map-ont $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_hybridassembly/out_spadesassembly/scaffolds.fasta $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_hybridassembly/herp_ont_trimmed.fastq.gz > reads.sam &&
    samtools view -S -b reads.sam | samtools sort - > reads.sorted.bam &&
    samtools index reads.sorted.bam


#ml BWA/0.7.17-intel-2018b

#bwa mem -t 10 -x ont2d canu_WGS_DH_hesp.contigs.fasta ../allreads.trimmed.fastq.gz > nanopore_mapping/canu_mapping.sam