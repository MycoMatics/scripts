cd #!/usr/bin/bash
#PBS -N nanopolish2.0.$PBS_JOBID
#PBS -l nodes=1:ppn=all
#PBS -o $PBS_JOBID.nanopolish2.0.stdout							
#PBS -e $PBS_JOBID.nanopolish2.0.stderr
#PBS -l walltime=48:00:00
#PBS -l mem=87gib
#PBS -m abe
#PBS -M pieter.asselman@ugent.be
    


# canu assembly
   #ml canu/2.1.1-GCCcore-10.2.0
   #cd $VSC_SCRATCH_VO_USER/minion/wgs_danny/assembly
    #canu -p canu_hespb -d canu genomeSize=17m -nanopore-raw $VSC_SCRATCH_VO_USER/minion/wgs_danny/trimmed/hespb_trimmed_q10_2000.fastq.gz useGrid=false

#############
#canu polish#
#############
  # Transform fastq to fasta
    #seqtk seq -a hespb_trimmed_q10_2000.fastq > herp_ontreads.fa
  # NANOPOLISH: indexing nanopore reads(needs fast5 files)
    ml parallel/20210322-GCCcore-10.2.0
    ml nanopolish/0.13.3-foss-2020b
    cd $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_polishassembly/nanopolish_2.0
      nanopolish index -d $VSC_SCRATCH_VO_USER/minion/wgs_danny/fast5/ -s $VSC_SCRATCH_VO_USER/minion/wgs_danny/fast5/sequencing_summary.txt ./herp_ontreads.fa &&
      
  # MINIMAP2-SAMTOOLS: index draft genome and align basecalled reads to draft
    #ml minimap2/2.24-GCCcore-11.2.0
    #ml SAMtools/1.14-GCC-11.2.0
    #cd $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_polishassembly/nanopolish_2.0
      #minimap2 -ax map-ont $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_polishassembly/nanopolish_2.0/canu_herp_15mb.contigs.fasta $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_polishassembly/nanopolish_2.0/herp_ontreads.fa > reads.sam
      #samtools view -S -b reads.sam | samtools sort - > reads.sorted.bam
      #samtools index reads.sorted.bam
  # NANOPOLISH: polish the draft genome
  # -P and -t needs to be adjusted so that P*t does not exceed the number of available threads
    #ml nanopolish/0.13.3-foss-2020b
    #ml parallel/20210322-GCCcore-10.2.0
    #cd $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_polishassembly/nanopolish_2.0
    nanopolish_makerange.py $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_polishassembly/nanopolish_2.0/canu_herp_15mb.contigs.fasta | parallel --results nanopolish.results -P 6 \
    nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_polishassembly/nanopolish_2.0/herp_ontreads.fa -b reads.sorted.bam -g $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_polishassembly/nanopolish_2.0/canu_herp_15mb.contigs.fasta -t 6 --min-candidate-frequency 0.1
    #nanopolish vcf2fasta -g $VSC_SCRATCH_VO_USER/minion/wgs_danny/project_polishassembly/nanopolish_2.0/canu_herp_15mb.contigs.fasta polished.*.vcf > polished_genome.fa
    
   # Note: this isn't actually an issue - it all works fine! I just thought I'd log what I did so if somebody searches for "slurm" in the Nanopolish issues, they'd find this.

#I use a cluster which is managed by Slurm and wanted to run Nanopolish as fast as possible. So instead of using parallel as shown in the Nanopolish README, I launched a separate job for each range:
#for range in $(python nanopolish_makerange.py draft.fa); do sbatch --nodes=1 --job-name=Nanopolish_"$range" --ntasks=1 --cpus-per-task=2 --mem=4096 --time=0-8:0:00 --wrap "nanopolish variants --consensus polished."$range.fa" -w "$range" -r reads.fa -b reads.sorted.bam -g draft.fa -t 2 --min-candidate-frequency 0.1"; done

#Then when they are all done, merge them together, as normal:
#python nanopolish_merge.py polished.*.fa > polished_genome.fa

#This let me run all 113 ranges of my bacterial genome simultaneously on different nodes of the cluster, and it was done in a couple of hours!