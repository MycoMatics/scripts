#!/usr/bin/bash
#PBS -N unicycler_timlep
#PBS -l nodes=1:ppn=half
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be


ml Unicycler/0.4.8-gompi-2020a-Python-3.8.2
cd $PBS_O_WORKDIR

unicycler -t 48 -s /scratch/gent/vo/001/gvo00142/vsc43352/timlep/illumina_raw/trimmomatic_out/Pailler206_R1_trimmed.fastq -l /scratch/gent/vo/001/gvo00142/vsc43352/timlep/dombea_0.6gb.fq.gz -o unicycler_hybrid_timlep


#unicycler -o $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_hybridassembly/out_unicyclersassembly \
#  -1 $VSC_SCRATCH_VO_USER/genomes/WMKJ_SRX7128323/trimmeddata/out_trimmomatic/SRR10432278_1_trimmed_paired.fastq \
# -2 $VSC_SCRATCH_VO_USER/genomes/WMKJ_SRX7128323/trimmeddata/out_trimmomatic/SRR10432278_2_trimmed_paired.fastq \
#  -l $VSC_SCRATCH_VO_USER/minion/wgs_danny/trimmed/hespb_trimmed_q10_2000.fastq.gz -t 96
  