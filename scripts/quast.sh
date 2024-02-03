#!/bin/bash
#SBATCH -J quast
#SBATCH --partition=long
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

Ref=$1
GFF=$2
ReadsF=$3
ReadsR=$4
ReadsONT=$5
OutDir=$6
Assembly=$7

apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/quast_5.2.0.sif quast.py \
    -t 16 -e -r $Ref -g $GFF \
    --k-mer-stats --circos \
    --pe1 $ReadsF --pe2 $ReadsR --nanopore $ReadsONT \
    -o $OutDir \
    $Assembly
