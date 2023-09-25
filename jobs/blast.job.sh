#!/usr/bin/bash

#PBS -N blast.$PBS_JOBID
#PBS -l nodes=1:ppn=8
#PBS -o $PBS_JOBID.blast.stdout							
#PBS -e $PBS_JOBID.blast.stderr
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be


ml BLAST+/2.14.0-gompi-2022b

cd $PBS_O_WORKDIR

for file in *.fa; do
 output_file="${file%.fa}.out" \
    blastn -db nt -query "$file" \
    -out "$output_file" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms qcovs" \
    -num_threads 20 \
done

#blastn -db /user/gent/433/vsc43352/scratch_vo/database/nt_ncbi/nt -query /user/gent/433/vsc43352/scratch_vo/minion/wgs_danny/project_contaminant_detect/herp_6000_filtlong_subsample.fasta -out ./results.sub_herp.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms qcovs" -num_threads 36
