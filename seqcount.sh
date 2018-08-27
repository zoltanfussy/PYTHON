#!/bin/sh
#PBS -N ryba-python
#PBS -l select=1:ncpus=1:mem=2gb:scratch_local=2gb
#PBS -l walltime=1:00:00
#PBS -m ae
#PBS -j oe

DATADIR="/storage/praha1/home/janko/Sterlet_transcriptome"

module add python27-modules-intel
cat $PBS_NODEFILE
cd $DATADIR

python seqcount.py