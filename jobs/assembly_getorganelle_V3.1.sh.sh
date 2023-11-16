#!/bin/bash

#PBS -N getorganelle
#PBS -l nodes=1:ppn=8
#PBS -o stdout.$PBS_JOBID						
#PBS -e stderr.$PBS_JOBID
#PBS -l walltime=01:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

#data
#ILLUMINA_TRIMMED=$VSC_SCRATCH_VO_USER/getorg/TLP/trimming/out_trimmomatic_11-36-16
ILLUMINA_TRIMMED=/data/gent/vo/001/gvo00142/vsc43352/cpgenomes/intro-cpgenomes/OX0001/trimmed-data
GETORG_OUT=/data/gent/vo/001/gvo00142/vsc43352/cpgenomes/intro-cpgenomes/getorganelle_out

# modules to load
ml GetOrganelle/1.7.5.3-foss-2021b


# create itteration trimmed_reads_sample file

for f in $ILLUMINA_TRIMMED/*_1_trimmed_paired.fastq
 do
 sample_name=$(basename "$f" _1_trimmed_paired.fastq)
 echo "$sample_name" >> $ILLUMINA_TRIMMED/samples_trimmed.txt #check if your files have the same extension layout! Adjust if necessary.
done &&

# load databases
get_organelle_config.py --add embplant_pt
#get_organelle_config.py -a all

#run getorg from reads
mkdir $GETORG_OUT/$files
while read files
do
	#for PE reads activate this line
	    get_organelle_from_reads.py -1 $ILLUMINA_TRIMMED/$files"_1_trimmed_paired.fastq" -2 $ILLUMINA_TRIMMED/$files"_2_trimmed_paired.fastq" -o $GETORG_OUT/$files --max-reads 3E7 -R 10 -k 21,45,65,85,105 -F embplant_pt -t 8
  #for SE reads activate this line
       #get_organelle_from_reads.py -u $ILLUMINA_TRIMMED/$files"_R1_trimmed.fastq" -o $GETORG_OUT/$files --max-reads 3E7 -R 25 -k 21,45,65,85,105 -F embplant_pt -t 40
#	echo "##########"
#	echo "PE sample_$files finished"
#	echo "##########"
done < $ILLUMINA_TRIMMED/samples_trimmed.txt