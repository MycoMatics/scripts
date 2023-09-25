#!/usr/bin/bash


#PBS -l walltime=05:00:00							
#PBS -N rpb2_allreads_nanoplot									
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l nodes=1:ppn=all							
#PBS -m abe
								


module load NanoComp/1.13.1-intel-2020b

cd $PBS_O_WORKDIR

for file in $PBS_O_WORKDIR/*.fq; do
  bn=`basename $file .fq`
  NanoComp -t 36 --fastq $file --names $bn -o ${bn}.nanocomp --plot violin
done

#NanoPlot -t 8 --fastq  <path to fastqfile> --loglength --N50 --title rpb2_minrun_summary -o ./nanoplot_rpb2

  # -t: number of threads used for this analysis, set at max avail in your system
  # --fastq: data format
  # --loglength: Additionally show logarithmic scaling of lengths in plots.
  # --N50 Show the N50 mark in the read length histogram
  # --title: Add a title to all plots, requires quoting if using spaces
  # -o: output folder
