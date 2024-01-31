#!/usr/bin/env bash
#SBATCH -J filtlong
#SBATCH --partition=medium
#SBATCH --mem=5G
#SBATCH --cpus-per-task=8

# porechop reads = $1

Reads=$1
Short=$(basename $1 _trim.fastq.gz)

filtlong --min_length 10000 --min_mean_q 93.7 $Reads | gzip > "$Short"_filt.fastq.gz

# percentage = 1 - (10^-(Q/10))
