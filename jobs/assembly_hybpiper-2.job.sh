{\rtf1\ansi\ansicpg1252\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 Times-Roman;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh14080\viewkind0
\deftab720
\pard\pardeftab720\sa240\partightenfactor0

\f0\fs24 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 #!/bin/bash\uc0\u8232 #PBS -N Hybpiper_assemble_BLASTx\u8232 #PBS -l nodes=1:ppn=63\u8232 #PBS -l walltime=72:00:00\
\'a0\
module load Hybpiper/2.3.2-foss-2024a\
\'a0\
unset OMP_PROC_BIND\
\'a0\
cd /scratch/gent/vo/001/gvo00142/vsc45815/Targeted_sequencing/Full_run\
\'a0\
reads_dir="trimmed_reads"\uc0\u8232 targets="/scratch/gent/vo/001/gvo00142/vsc45815/Targeted_sequencing/Target_sequences/target_file_milkcap291_full_genes_Russulales_AA.fasta"\
\'a0\
mkdir -p /tmp/hybpiper_db\uc0\u8232 cp "$targets" /tmp/hybpiper_db/targets_local.fasta\u8232 local_targets="/tmp/hybpiper_db/targets_local.fasta"\u8232 echo "Using local target file: $local_targets"\u8232 \u8232 for f in $\{reads_dir\}/*_1_paired.fq.gz; do\u8232 \'a0\'a0\'a0 sample=$(basename "$f" _1_paired.fq.gz)\u8232 \'a0\'a0\'a0 echo "Running HybPiper for sample: $sample"\u8232 \'a0\'a0\'a0 hybpiper assemble \\\u8232 \'a0\'a0\'a0\'a0\'a0\'a0\'a0 -r $\{reads_dir\}/$\{sample\}_1_paired.fq.gz $\{reads_dir\}/$\{sample\}_2_paired.fq.gz \\\u8232 \'a0\'a0\'a0\'a0\'a0\'a0\'a0 -t_aa $local_targets \\\u8232 \'a0\'a0\'a0\'a0\'a0\'a0\'a0 --prefix $sample \\\u8232 \u8195 \u8195 \u8195 \u8195 --compress_sample_folder\
}