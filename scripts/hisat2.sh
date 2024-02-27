#!/usr/bin/env bash
#SBATCH -J hisat2
#SBATCH --partition=medium
#SBATCH --mem=18G
#SBATCH --cpus-per-task=4

indexDir=$1
F_Reads=$2
R_Reads=$3

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate rna-seq

Short=$(basename $2 .R1_val_1.fq.gz.sam)

hisat2 --dta -p 8 -x $indexDir -1 $F_Reads -2 $R_Reads -S $Short.sam --summary-file $Short.summary

samtools sort -@ 8 -o $Short.sorted.bam $Short.sam

rm $Short.sam
