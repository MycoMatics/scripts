#!/bin/bash

#PBS -N $PBS_JOBID.TLP_getorganelle
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

#data
#ILLUMINA_SE_TRIMMED=$VSC_SCRATCH_VO_USER/getorg/TLP/trimming/out_trimmomatic_11-36-16
ILLUMINA_PE_TRIMMED=$VSC_SCRATCH_VO_USER/getorg/test_script/PE_illumina
GETORG_OUT=$VSC_SCRATCH_VO_USER/getorg/TLP/getorg_output

# modules to load
module load GetOrganelle/1.7.4-pre2-foss-2020b

# create itteration trimmed_reads_sample file

for f in $ILLUMINA_SE_TRIMMED/*.fastq
 do
  echo "$(basename "$f" "_R1_trimmed.fastq")" #check if your files have the same extension layout! Adjust if necessary.
done > $ILLUMINA_SE_TRIMMED/samples_trimmed.txt &&

# load databases
get_organelle_config.py --add embplant_pt

#run getorg from reads

while read files
do
	#for PE reads activate this line
   mkdir $GETORG_OUT/$files
	    get_organelle_from_reads.py -1 $ILLUMINA_PE_TRIMMED/$files".1.fq.gz" -2 $ILLUMINA_PE_TRIMMED/$files".2.fq.gz" -o $GETORG_OUT/$files --max-reads 3E7 -R 10 -k 21,45,65,85,105 -F embplant_pt -t 8
  #for SE reads activate this line
       #get_organelle_from_reads.py -u $ILLUMINA_SE_TRIMMED/$files"_R1_trimmed.fastq" -o $GETORG_OUT/$files --max-reads 3E7 -R 25 -k 21,45,65,85,105 -F embplant_pt -t 40
	echo "##########"
	echo "PE sample_$files finished"
	echo "##########"
done < $ILLUMINA_SE_TRIMMED/samples_trimmed.txt
