
#!/bin/bash
#PBS -N cutadapt.$PBS_JOBID
#PBS -l nodes=1:ppn=8
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=0:30:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

ml cutadapt/3.5-GCCcore-11.2.0

 
Cutadapt
				$ cutadapt -j 8 -q 30 -m 50 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o reads_1_trimmed.fastq -p reads_2_trimmed.fastq rawreads_1.fastq.gz rawreads_2.fastq.gz

					-m LENGTH: Discard trimmed reads that are shorter than LENGTH.
					-a 3'/3 adapter to be removed from the first read in a pair.
					-A 3'/3 adapter to be removed the second read in a pair.
					-o FILE: Write trimmed first reads to FILE.
					-p FILE: Write trimmed second reads to FILE.
					-j CORES: Number of CPU cores to use.