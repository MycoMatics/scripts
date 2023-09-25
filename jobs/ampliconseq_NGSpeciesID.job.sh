#!/bin/bash

#PBS -l walltime=48:00:00							
#PBS -N rpb2_NGSpeciesID.$PBS_JOBID									
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l nodes=1:ppn=all							
#PBS -m abe
								


module load NGSpeciesID/0.1.2.1-foss-2021b

cd $PBS_O_WORKDIR

for file in <path to>/*.fastq; do 
  bn=`basename $file .fastq`
  NGSpeciesID --ont --t 36 --consensus --medaka --fastq $file --outfolder ${bn}
done



# --ont: reads in fastq format are from Oxford nanopore reads
# --t: number of threads to use, set according to available ppn on your system
# --consensus: We need a consensussequence per sample
# --racon: all reads will be polished (mapped to the consensus) by the tool racon
# --racon_iter: iterations of polishing
# --fastq: input file
# --outfolder output folder per sample

#NGSpeciesID --ont --m 650 --s 100 --fastq sample_7552.fastq --outfolder ./run1 --t 10 --consensus --medaka
#NGSpeciesID --ont --m 650 --s 100 --fastq sample_7552.fastq --outfolder ./run2 --t 10 --consensus --racon --racon_iter 5
#--abundance_ratio 0,01 add this flag to build consensus from low abundant reads
#--m TARGET_LENGTH     Intended amplicon length. Invoked to filter out reads with length greater than m + s or smaller than m - s (default = 0 which means no filtering) (default: 0)
#--s TARGET_DEVIATION  Maximum allowed amplicon-length deviation. Invoked to filter out reads with length greater than m + s or smaller than m - s (default = 0 which means nofiltering) (default: 0)