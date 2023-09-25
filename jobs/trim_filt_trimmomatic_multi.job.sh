
#!/bin/bash
#PBS -N trimmomatic.$PBS_JOBID
#PBS -l nodes=1:ppn=8
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=0:30:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be


# Data Paths

ILLUMINA_RAWDATA=/user/gent/433/vsc43352/scratch_vo/genomes/WMKJ_SRX7128323/rawdata #take notice of the file extension names different options '.fq.gz' '.fastq.gz'
ILLUMINA_ADAPERS=/user/gent/433/vsc43352/scripts

# output directories qc
#TRIMMOMATIC_OUT=/user/gent/433/vsc43352/scratch_vo/genomes/WMKJ_SRX7128323/rawdata/out_trimmomatic_$(date +"%Y-%m-%d_%H-%M-%S")
TRIMMOMATIC_OUT=/user/gent/433/vsc43352/scratch_vo/genomes/WMKJ_SRX7128323/rawdata/out_trimmomatic


# Environmentals
#VSC_DATA=/data/gent/433/vsc43352
#VSC_DATA_VO=/data/gent/vo/001/gvo00142
#VSC_DATA_VO_GENT=/data/gent/vo/001/gvo00142
#VSC_DATA_VO_USER=/data/gent/vo/001/gvo00142/vsc43352
#VSC_DATA_VO_USER_GENT=/data/gent/vo/001/gvo00142/vsc43352
#VSC_HOME=/user/gent/433/vsc43352
#VSC_SCRATCH_ARCANINE=/arcanine/scratch/gent/433/vsc43352
#VSC_SCRATCH_ARCANINE_VO=/arcanine/scratch/gent/vo/001/gvo00142
#VSC_SCRATCH_ARCANINE_VO_USER=/arcanine/scratch/gent/vo/001/gvo00142/vsc43352
#VSC_SCRATCH_KYUKON=/kyukon/scratch/gent/433/vsc43352
#VSC_SCRATCH_KYUKON_VO=/kyukon/scratch/gent/vo/001/gvo00142
#VSC_SCRATCH_KYUKON_VO_USER=/kyukon/scratch/gent/vo/001/gvo00142/vsc4335
#VSC_SCRATCH_VO=/scratch/gent/vo/001/gvo00142
#VSC_SCRATCH_VO_USER=/scratch/gent/vo/001/gvo00142/vsc43352
#VSC_VO_GENT=gvo00142
#VSC_VO=gvo00142


# QC directories
mkdir $TRIMMOMATIC_OUT

# Create sampleslist to itterate, serves as input for itteration process in trimmomatic
for f in $ILLUMINA_RAWDATA/*_1.fastq.gz
 do
  echo "$(basename "$f" "_1.fastq.gz")" #check if your files have the same extension layout! Adjust if necessary.
done > $ILLUMINA_RAWDATA/samples.txt &&


# Load modules
module load Trimmomatic/0.39-Java-11

### Run trimmomatic SE OR PE reads #-out what you don't need!

while read p
do 

### ACTIVATE FOR SE data
#java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -phred33 -trimlog $TRIMMOMATIC_OUT/$p".log" $ILLUMINA_RAWDATA/$p"_R1.fastq.gz" $TRIMMOMATIC_OUT/$p"_R1_trimmed.fastq" ILLUMINACLIP:$ILLUMINA_ADAPERS/alladapterstrimmomatic.fa:2:30:10:1:TRUE SLIDINGWINDOW:5:20

### ACTIVATE FOR PE data
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -trimlog $TRIMMOMATIC_OUT/$p".log" \
      $ILLUMINA_RAWDATA/$p"_1.fastq.gz" \
      $ILLUMINA_RAWDATA/$p"_2.fastq.gz" \
      $TRIMMOMATIC_OUT/$p"_1_trimmed_paired.fastq" \
      $TRIMMOMATIC_OUT/$p"_1_trimmed_unpaired.fastq" \
      $TRIMMOMATIC_OUT/$p"_2_trimmed_paired.fastq" \
      $TRIMMOMATIC_OUT/$p"_2_trimmed_unpaired.fastq" \
      ILLUMINACLIP:$ILLUMINA_ADAPERS/alladapterstrimmomatic.fa:2:30:10:1:TRUE SLIDINGWINDOW:5:20
  
done < $ILLUMINA_RAWDATA/samples.txt