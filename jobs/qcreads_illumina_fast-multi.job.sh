#!/bin/bash

#PBS -N $PBS_JOBID.fastqc
#PBS -l nodes=1:ppn=8
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=00:30:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be


# raw datapath
ILLUMINA=//user/gent/433/vsc43352/scratch_vo/genomes/WMKJ_SRX7128323/rawdata/out_trimmomatic

# output directories qc
FASTQC_OUT=//user/gent/433/vsc43352/scratch_vo/genomes/WMKJ_SRX7128323/rawdata/out_fastqc_trimmed
MULTIQC_OUT=//user/gent/433/vsc43352/scratch_vo/genomes/WMKJ_SRX7128323/rawdata/out_multiqc_trimmed

#environmentals
#VSC_DATA=/data/gent/433/vsc43352
#VSC_DATA_VO=/data/gent/vo/001/gvo00142
#VSC_DATA_VO_GENT=/data/gent/vo/001/gvo00142
#VSC_DATA_VO_USER=/data/gent/vo/001/gvo00142/vsc43352
#VSC_DATA_VO_USER_GENT=/data/gent/vo/001/gvo00142/vsc43352
#VSC_HOME=/user/gent/433/vsc43352
#VSC_SCRATCH_ARCANINE=/arcanine/scratch/gent/433/vsc43352
#VSC_SCRATCH_ARCANINE_VO=/arcanine/scratch/gent/vo/001/gvo00142
#VSC_SCRATCH_ARCANINE_VO_USER=/arcanine/scratch/gent/vo/001/gvo00142/vsc43352
#VSC_SCRATCH_KYUKON=/kyukon/scratch/gent/433/vsc43352
#VSC_SCRATCH_KYUKON_VO=/kyukon/scratch/gent/vo/001/gvo00142
#VSC_SCRATCH_KYUKON_VO_USER=/kyukon/scratch/gent/vo/001/gvo00142/vsc4335
#VSC_SCRATCH_VO=/scratch/gent/vo/001/gvo00142
#VSC_SCRATCH_VO_USER=/scratch/gent/vo/001/gvo00142/vsc43352
#VSC_VO_GENT=gvo00142
#VSC_VO=gvo00142


# QC directories
mkdir $FASTQC_OUT
mkdir $MULTIQC_OUT

### run fastQC
module load FastQC/0.11.9-Java-11

fastqc $ILLUMINA/*.fastq  -o $FASTQC_OUT && module del FastQC/0.11.9-Java-11 && module load MultiQC/1.9-intel-2020a-Python-3.8.2 && multiqc $FASTQC_OUT -o $MULTIQC_OUT 
