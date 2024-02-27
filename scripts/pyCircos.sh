#!/usr/bin/env bash
#SBATCH -J circos
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

python /mnt/shared/scratch/jnprice/varnes_genome/scripts/pyCircos_test.py