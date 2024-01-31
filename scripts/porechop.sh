#!/usr/bin/env bash
#SBATCH -J porechop
#SBATCH --partition=himem
#SBATCH --mem=300G
#SBATCH --cpus-per-task=4

# reads_dir = $1
# output = $2

porechop \
    -t 8 \
    -i $1 \
    -o $2
