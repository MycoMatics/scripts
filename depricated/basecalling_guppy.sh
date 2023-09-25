#!/bin/bash

#PBS -S /bin/bash
#PBS -l walltime=06:00:00							
#PBS -N guppy.$PBS_JOBID									
#PBS -o $PBS_JOBID.guppy.stdout							
#PBS -e $PBS_JOBID.guppy.stderr
#PBS -l nodes=1:ppn=all:gpus=2							
#PBS -m abe
#PBS -M pieter.asselman@ugent.be								


# DATArun on 2GPU/ 11Gbases - ~3h
 
## BASECALL ##
# guppy_basecaller -i <input FAST5 FOLDER> -s <OUTPUTFOLDER> -c dna_r10.3_450bps_hac.cfg --compress_fastq --verbose_logs -q 0 -x 'cuda:all:100%'
	#-i 			Path to input fast5 files
	#-s 			Path to save fastq files
	# -c 			config file to use
		## dna_r10.3_450bps_sup.cfg
		## /opt/ont/guppy/data (holds jsn and cfg files)
	# --compres_fastq 	Compress fastq output files with gzip
	# --verbose_logs	Enable verbose logs
	# -q 0			Maximum number of records per fastq file, 0 means use a single file (per worker, per run id)
	# -x 'cuda:all:100%	Specify basecalling device, all 100% means all GPU resources

## DEMULTIPLEX only ##
# guppy_barcoder -i <basecalled fastq files> -s <folder to demultiplexed data> --barcode_kits EXP-PBC001 --compress_fastq --trim_barcodes --num_extra_bases_trim -3 --verbose_logs -q 0 -x 'cuda:all:100%'

## BASECALL and DEMULTIPLEX ##
# guppy_basecaller -i <input FAST5 FOLDER> -s <folder to demultiplexed data> -c dna_r10.3_450bps_hac.cfg --barcode_kits EXP-PBC001 --compress_fastq --trim_barcodes --num_extra_bases_trim -3 --verbose_logs -q 0 -x 'cuda:all:100%'


module load Guppy/4.4.1-gpu
module load CUDA/11.0.2-GCC-9.3.0


guppy_barcoder -i /scratch/gent/vo/000/gvo00058/vsc43352/minion_runs/20210222_ITSBC1_12/2102022_ITSBC1_12/20210222_1452_MN35631_FAP52007_80ea3e1d/fastq -s /data/gent/vo/000/gvo00058/vsc43352/minion_runs/20210222_ITSBC1_12/fastq/demultiplexed --barcode_kits EXP-PBC001 --compress_fastq --trim_barcodes --num_extra_bases_trim -3 --verbose_logs -q 0 -x 'cuda:all:100%'

