#!/usr/bin/bash


#PBS -l walltime=05:00:00							
#PBS -N nanofilt									
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l nodes=1:ppn=all							
#PBS -m abe
								


module load nanofilt/2.6.0-intel-2020a-Python-3.8.2

cd $PBS_O_WORKDIR

cat <path to your fasq file> | NanoFilt -q 12 -l 800 --maxlength 1040 > filtered_reads.fq

#ITTERATE
	# for file in *.fastq; do
	# bn=`basename $file .fastq`
	# cat $file | NanoFilt -q 12 -l 400 --maxlength 1300 > ../../4ngspeciesid/${bn}_filtered.fastq
  #done
	
#EXPLANATION
  # cat: regex for view content of a file
  # | : contente of file will be written to stdout via the terminal, the | takes the output to the next command
  # -q: cut off value for min quality score (q10 ~ 90% accuracy
              #q10 - 1 in 100 - 90%
              #q20 - 1 in 100 - 99%
              #q30 - 1 in 1000 - 99.9%
  # -l: cut off value for min length
  # --maxlength: cut off value for max length
  # >: redirect output from NanoFilt to a file
              