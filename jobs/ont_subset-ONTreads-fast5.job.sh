#!/usr/bin/bash
#PBS -N slaking_interactive.$PBS_JOBID
#PBS -l nodes=1:ppn=8
#PBS -l mem=26GB
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=05:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

module load ont-fast5-api/3.3.0-foss-2020b


##go to jobfolder, this is the flolder from where you call and launch your script

cd $PBS_O_WORKDIR


fast5_subset --input $VSC_DATA_VO_USER/HPC-ampseq-course/dataset/fast5/fast5 --save_path $VSC_DATA_VO_USER/HPC-ampseq-course/dataset/fast5_subset --read_id_list sequencing_summary_barcode01_only_andheaders.txt --batch_size 1000 --recursive