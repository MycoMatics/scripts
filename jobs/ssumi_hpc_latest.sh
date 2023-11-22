#!/bin/bash

  # "name" of the job (optional)
#PBS -N ssUMI_pipeline_all_barcodes

  # specification (required!)
  #   nodes=   number of nodes; 1 for serial; 1 or more for parallel
  #   ppn=     number of processors per node; 1 for serial; up to 8
  #   if you want your "private" node: ppn=8
  #   mem=     memory required
  
  #per UMI barcode of 1GB > 1hour is needed (usearch can only use 10 threads)
#PBS -l walltime=40:00:00
#PBS -l nodes=1:ppn=all

  # send mail notification (optional)
  #   a        when job is aborted
  #   b        when job begins
  #   e        when job ends
  #   M        your e-mail address (should always be specified)

#PBS -m a
#PBS -m b
#PBS -m e
#PBS -M glen.dierickx@ugent.be

  # redirect standard output (-o) and error (-e) (optional)
  # if omitted, the name of the job (specified by -N) or
  # a generic name (name of the script followed by .o or .e and 
  # job number) will be used
#PBS -o stdout.$PBS_JOBID
#PBS -e stderr.$PBS_JOBID

# go to the (current) working directory (optional, if this is the directory where you submitted the job)
#cd path to files or
cd $PBS_O_WORKDIR

#execute the script
echo "Start job - Today is $(date +%Y-%m-%d)"
echo "amount of virtual memory Kb that is available: $(ulimit -v)"
#load modules, activate environments
source activate
conda activate longread_umi
echo "this is conda environment $CONDA_DEFAULT_ENV"
#set variables
FORWARD_ADAPTER="GTATCGTGTAGAGACTGCGTAGG" #-f flag, nanopore adapter
FORWARD_PRIMER="TCCGTAGGTGAACCTGCGG" #-F flag,ITS1
REVERSE_ADAPTER="AGTGATCGAGTCAGTGCGAGTG" #-r flag, nanopore adapter
REVERSE_PRIMER="TCCTCCGCTTATTGATATGC" #-R flag, ITS4
MODEL="r1041_e82_400bps_sup_v4.2.0" #specify basecalling model used
THREADS="20" # -T flag, set depending on submitted cluster
RACON_ROUNDS="3" #-c FLAG
MEDAKA_ROUNDS="2" #-p flag
MEDAKA_THREADS=$(expr $THREADS / 5) #-T flag

#run script
for FILE in *.fastq; do
	echo starting with $(basename "$FILE")
	longread_umi ssumi_std \
	  -d $FILE \
	  -v 3 \
	  -o ${FILE%.fastq}_ssumi_pipe \
	  -s 200 \
	  -e 200 \
	  -E 0.2 \
	  -m 450 \
	  -M 1200 \
	  -f $FORWARD_ADAPTER \
	  -F $FORWARD_PRIMER \
	  -r $REVERSE_ADAPTER \
	  -R $REVERSE_PRIMER \
	  -c $RACON_ROUNDS \
	  -p $MEDAKA_ROUNDS \
	  -q $MODEL \
	  -t $THREADS \
	  -T $MEDAKA_THREADS
	echo done with UMIs of $(basename "$FILE")
done
#clean-up
module purge
echo \############ pipeline done \##############

