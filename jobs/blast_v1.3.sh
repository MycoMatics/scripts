#!/usr/bin/bash

#PBS -N blast.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o $PBS_JOBID.blast.stdout							
#PBS -e $PBS_JOBID.blast.stderr
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be


ml BLAST+/2.14.0-gompi-2022b

cd $PBS_O_WORKDIR

# Available Data bases installed. => !!! ADD PATH TO # VARIABLES !!!

    # NCBI NT database: /data/gent/vo/001/gvo00142/data_share_group/databases_blast/nt/nt
    # UNITE database: /data/gent/vo/001/gvo00142/data_share_group/databases_blast/unite/unite


# VARIABLES
DB=/data/gent/vo/001/gvo00142/data_share_group/databases_blast/unite/unite
IDENTITY=90
MAX_TARGET_SEQS=1
THREADS=20
OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms qcovs"

for file in *.fa; do
 output_file="${file%.fa}.blast.out"

  blastn -db "$DB" \
 -query "$file" \
 -out "$output_file" \
 -outfmt "$OUTFMT" \
 -perc_identity "$IDENTITY" \
 -max_target_seqs "$MAX_TARGET_SEQS" \
 -num_threads "$THREADS"
done