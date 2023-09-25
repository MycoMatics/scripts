
#!/bin/bash
#PBS -N trimmomatic.$PBS_JOBID
#PBS -l nodes=1:ppn=8
#PBS -l mem=26GiB
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=0:30:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be


module load Trimmomatic/0.39-Java-11

cd $PBS_O_WORKDIR

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
	210601_I23_V350004411_L4_OX0001_1.fq.gz \
	210601_I23_V350004411_L4_OX0001_1.fq.gz \
	OX1_1_paired.fq.gz \
	OX1_1_unpaired.fq.gz \
	OX1_2_paired.fq.gz \
	OX1_2_unpaired.fq.gz \
	ILLUMINACLIP:alladapterstrimmomatic.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15 
 
