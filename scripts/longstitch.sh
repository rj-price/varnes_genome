#!/bin/bash
#SBATCH -J longstitch
#SBATCH --partition=medium
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8

#! Files need to be in same directory
#! Assembly needs to have .fa extension
#! Reads need to have .fq.gz extension

# Assembly = $1
# ONT Reads = $2
# Genome Size = $3

assemblyShort=$(basename $1 .fasta)
readsShort=$(basename $2 .fastq.gz)
genomeSize=$3

cp $1 ./$assemblyShort.fa
cp $2 ./$readsShort.fq.gz

apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/longstitch_1.0.5.sif \
    longstitch run t=16 draft=$assemblyShort reads=$readsShort G=$genomeSize outprefix=longstitch
