#!/usr/bin/bash
#PBS -N minibar_rpb2.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be


module load minibar/20200326-iccifort-2020.1.217-Python-3.8.2

cd $PBS_O_WORKDIR
minibar.py <path to barcode file>.txt <pah to filtered reads>.fastq -F -T -e 5 -E 8


#ITTERATE
	#file= barcodefile: IndexCombinationMyco.txt
	#file2= reads: allreads_filtered.fastq

	#for file in *.txt; do
	#file2="${file%.txt}_allreads_filtered.fastq
	#minibar.py $file $file2 -F -T -e 5 -E 8
	#done

#ONELINE CODE
	#for file in *.txt; do file2="${file%.txt}_allreads_filtered.fastq"; minibar.py $file $file2 -F -T -e 5 -E 8; done


#EXPLANATION
  # -e: Index edit dist 5
  # -E: Primer edit dist 8
  # -F: create individual sample files for sequences
  # -T: output type, trims barcode and primer from each end of the sequence, then outputs record

