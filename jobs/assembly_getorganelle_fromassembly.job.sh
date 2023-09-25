#!/usr/bin/bash
#PBS -N getorg_fromONTassembly
#PBS -l nodes=1:ppn=half
#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID
#PBS -o stderr.$PBS_JOBNAME.$PBS_JOBID
#PBS -l walltime=24:00:00
#PBS -m abe
#PBS -M pieter.asselman@ugent.be

module load GetOrganelle/1.7.5.3-foss-2021b 
ml Bandage/0.9.0-GCCcore-11.2.0

cd $PBS_O_WORKDIR

get_organelle_from_assembly.py -F embplant_pt --min-depth 50 -g assembly_graph_flye.gfa -o getorg_flye_out -t 48 --verbose


