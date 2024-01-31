#!/usr/bin/env bash
#SBATCH -J filtlong
#SBATCH --partition=medium
#SBATCH --mem=5G
#SBATCH --cpus-per-task=8

# porechop reads = $1

Reads=$1
Short=$(basename $1 _trim.fastq.gz)

filtlong --min_length 1000 --min_mean_q 98 $Reads | gzip > "$Short"_highqual.fastq.gz

# Q17
# percentage = 1 - (10^-(Q/10))
