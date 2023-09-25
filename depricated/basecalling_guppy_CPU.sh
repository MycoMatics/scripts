#!/bin/bash

#PBS -l walltime=01:00:00							
#PBS -N guppy.$PBS_JOBID									
#PBS -o guppy.stdout.$PBS_JOBID							
#PBS -e guppy.stderr.$PBS_JOBID
#PBS -l nodes=1:ppn=all							
#PBS -m abe
#PBS -M pieter.asselman@ugent.be								

#ml Guppy/6.1.7-gpu
ml Guppy/6.1.7-cpu
ml CUDA/11.4.1

cd $PBS_O_WORKDIR


# DATArun on 2GPU/ 11Gbases - ~3h
 
## BASECALL ##
guppy_basecaller -i $PBS_O_WORKDIR -s $PBS_O_WORKDIR/CPU_basecalling.$PBS_JOBID -c dna_r9.4.1_450bps_fast.cfg --compress_fastq --verbose_logs -q 0 --num_callers 36
	#-i 			Path to input fast5 files
	#-s 			Path to save fastq files
	# -c 			config file to use
		## dna_r10.3_450bps_sup.cfg
		## /opt/ont/guppy/data (holds jsn and cfg files)
	# --compres_fastq 	Compress fastq output files with gzip
	# --verbose_logs	Enable verbose logs
	# -q 0			Maximum number of records per fastq file, 0 means use a single file (per worker, per run id)
	# -x 'cuda:all:100%	Specify basecalling device, all 100% means all GPU resources
	# --chunks_per_caller: A soft limit on the number of chunks in each basecaller's chunk queue. When a read is sent to the basecaller, it is broken up into “chunks” of signal, and each chunk is basecalled in isolation. Once                           all the chunks for a read have been basecalled, they are combined to produce a full basecall. --chunks_per_caller sets a limit on how many chunks will be collected before they are dispatched for                                 basecalling. On GPU platforms this is an important parameter to obtain good performance, as it directly influences how much computation can be done in parallel by a single basecaller.
  # Number of parallel callers (--num_callers): Number of parallel basecallers to create. A thread will be spawned for each basecaller to use. Increasing this number will allow Guppy to make better use of multi-core CPU                                                       systems, but may impact overall system performance.

## DEMULTIPLEX only ##
# guppy_barcoder -i <basecalled fastq files> -s <folder to demultiplexed data> --barcode_kits EXP-PBC001 --compress_fastq --trim_barcodes --num_extra_bases_trim -3 --verbose_logs -q 0 -x 'cuda:all:100%'

## BASECALL and DEMULTIPLEX ##
# guppy_basecaller -i <input FAST5 FOLDER> -s <folder to demultiplexed data> -c dna_r10.3_450bps_hac.cfg --barcode_kits EXP-PBC001 --compress_fastq --trim_barcodes --num_extra_bases_trim -3 --verbose_logs -q 0 -x 'cuda:all:100%'


#module load Guppy/4.4.1-gpu
#module load CUDA/11.0.2-GCC-9.3.0


#guppy_barcoder -i /scratch/gent/vo/000/gvo00058/vsc43352/minion_runs/20210222_ITSBC1_12/2102022_ITSBC1_12/20210222_1452_MN35631_FAP52007_80ea3e1d/fastq -s /data/gent/vo/000/gvo00058/vsc43352/minion_runs/20210222_ITSBC1_12/fastq/demultiplexed --barcode_kits EXP-PBC001 --compress_fastq --trim_barcodes --num_extra_bases_trim -3 --verbose_logs -q 0 -x 'cuda:all:100%'

