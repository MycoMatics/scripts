#!/usr/bin/bash
#PBS -N hybpiper_conda_tutorialdata
#PBS -l nodes=1:ppn=8
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=01:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

cd $PBS_O_WORKDIR

source ~/.bashrc
conda activate hybpiper &&

# Unpack the test dataset
tar -zxf test_reads.fastq.tar.gz

# Remove any previous runs
parallel rm -r {} :::: namelist.txt


# Run main HybPiper command with all available CPUs
while read sample_name
do
  hybpiper assemble -r ${sample_name}*.fastq -t_dna test_targets.fasta --prefix ${sample_name} --bwa  --run_intronerate
done < namelist.txt


# Get runs statistics
hybpiper stats -t_dna test_targets.fasta gene namelist.txt


# Get heatmap of length recovery
hybpiper recovery_heatmap seq_lengths.tsv

# Recover DNA and amino-acid sequences
hybpiper retrieve_sequences -t_dna test_targets.fasta dna --sample_names namelist.txt --fasta_dir 01_dna_seqs
hybpiper retrieve_sequences -t_dna test_targets.fasta aa --sample_names namelist.txt --fasta_dir 02_aa_seqs


# Recover paralog sequences
hybpiper paralog_retriever namelist.txt -t_dna test_targets.fasta


echo "DONE!"

