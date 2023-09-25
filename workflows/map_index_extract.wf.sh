#!/usr/bin/bash
#PBS -N map2ref_hespb_WMKJ01.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -e stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=48:00:00
#PBS -l mem=100gb
#PBS -m abe
#PBS -M pieter.asselman@ugent.be





ml minimap2/2.20-GCCcore-10.2.0
ml SAMtools/1.11-GCC-10.2.0
ml Qualimap/2.2.1-foss-2020b-R-4.0.3
ml KAT/2.4.2-foss-2019a-Python-3.7.2

cd $VSC_SCRATCH_VO_USER/minion/wgs_danny/map2ref

# SHORT(map,bam,index)
 	#minimap2 -ax map-ont -t 8 ref.fa reads.fasta | samtools sort -o reads-ref.sorted.bam -T reads.tmp
 	#samtools index reads-ref.sorted.bam
# minimap - map to ref
	#echo 'START minimap - map to ref' &&
	#minimap2 -ax map-ont WMKJ01.1.DH.fasta.gz ../trimmed/hespb_trimmed_q10_2000.fastq.gz > hespb_WMKJ01_mapped.sam
	#echo 'DONE map to ref' &&

# sam to bam
	#echo 'START sam to bam'
	#samtools view -S -b hespb_WMKJ01_mapped.sam | samtools sort - > hespb_WMKJ01_mapped_sort.bam &&
	#echo 'DONE sam to bam'

# index bam
	#echo 'START index bam' &&
	#samtools index hespb_WMKJ01_mapped_sort.bam
	#echo 'DONE index bam' &&

# stats
	#echo 'START statistics mapping flagstat'
	#samtools flagstat hespb_WMKJ01_mapped_sort.bam &&
	#echo 'DONE statistics mapping flagstat'

# qualimap
	echo 'START qualimap' &&
	qualimap bamqc -bam ./hespb_WMKJ01_mapped_sort.bam -outdir hesp_WMKJ01_mapped_sort_bamqc 
	echo 'DONE qualimap' &&

# extract mapped reads from bam
	echo 'START extracting mapped reads from bam'	
	samtools view -b -F 4 hespb_WMKJ01_mapped_sort.bam > hespb_WMKJ01_mappedreadsonly.bam &&
	echo 'DONE mapped reads from bam extracted'

# bam to fastq
	echo 'START convert bam to fastq' &&
	samtools fastq hespb_WMKJ01_mappedreadsonly.bam > hespb_WMKJ01_mappedreadsonly.fastq
	echo 'DONE mapped read only fastq generated' &&

