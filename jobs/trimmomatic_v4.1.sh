#!/bin/bash
#PBS -N trimmomatic
#PBS -l nodes=1:ppn=8
#PBS -o stdout.$PBS_JOBID						
#PBS -e stderr.$PBS_JOBID
#PBS -l walltime=0:30:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

# Data Paths
ILLUMINA_RAWDATA=/yourpathto/OX0001 #take notice of the file extension names different options '.fq.gz' '.fastq.gz'
ILLUMINA_ADAPTERS=/yourpathto/adapters

# make trimmed data directory
mkdir $ILLUMINA_RAWDATA/trimmed-data

# Create sampleslist to itterate, serves as input for itteration process in trimmomatic
for file  in $ILLUMINA_RAWDATA/*_1.fq.gz
 do
  sample_name=$(basename "$file" _1.fq.gz) # extract sample name
  echo "${sample_name}" >> $ILLUMINA_RAWDATA/samples.txt
  done &&


# Load modules
module load Trimmomatic/0.39-Java-11

### Run trimmomatic SE OR PE reads #-out what you don't need!

while read p
do 

### ACTIVATE FOR SE data
#java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -phred33 \
#-trimlog $TRIMMOMATIC_OUT/$p".log" \
#$ILLUMINA_RAWDATA/$p"_R1.fastq.gz" \
#$TRIMMOMATIC_OUT/$p"_R1_trimmed.fastq" \
#ILLUMINACLIP:$ILLUMINA_ADAPERS/alladapterstrimmomatic.fa:2:30:10:1:TRUE SLIDINGWINDOW:5:20

### ACTIVATE FOR PE data
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -trimlog $ILLUMINA_RAWDATA/trimmed-data/$p".log" \
      $ILLUMINA_RAWDATA/$p"_1.fq.gz" \
      $ILLUMINA_RAWDATA/$p"_2.fq.gz" \
      $ILLUMINA_RAWDATA/trimmed-data/$p"_1_trimmed_paired.fastq" \
      $ILLUMINA_RAWDATA/trimmed-data/$p"_1_trimmed_unpaired.fastq" \
      $ILLUMINA_RAWDATA/trimmed-data/$p"_2_trimmed_paired.fastq" \
      $ILLUMINA_RAWDATA/trimmed-data/$p"_2_trimmed_unpaired.fastq" \
      ILLUMINACLIP:$ILLUMINA_ADAPTERS/alladapterstrimmomatic.fa:2:30:10:1:TRUE SLIDINGWINDOW:5:20
  
done < $ILLUMINA_RAWDATA/samples.txt
