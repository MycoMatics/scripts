#!/bin/bash
#!/bin/bash
#PBS -N trimmomatic.$PBS_JOBID
#PBS -l nodes=1:ppn=8
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=0:30:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

# Load modules
module load Trimmomatic/0.39-Java-11

while read sample
do
	echo "#######"
	echo "$sample"
	trimmomatic PE -phred33 \
	-trimlog trimmed/$ample'.log' \
	pathto/$sample'_R1.fastq.gz' pathto/$sample'_R2.fastq.gz' \
	destinationpath/$sample'_R1_paired.fastq' destinationpath/$sample'_R1_unpaired.fastq' destinationpath/$sample'_R2_paired.fastq'\
	destinationpath/$sample'_R2_unpaired.fastq' ILLUMINACLIP:/pathto/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE SLIDINGWINDOW:5:20 \
done < samples.txt	

