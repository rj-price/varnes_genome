#!/usr/bin/env bash
#SBATCH -J satsumasynteny
#SBATCH --partition=medium
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

# INPUT
targetGenome=$1
queryGenome=$2
outDir=$3

export SATSUMA2_PATH=/mnt/shared/scratch/jnprice/apps/conda/pkgs/satsuma2-20161123-h9f5acd7_4/bin

SatsumaSynteny2 -threads 8 -km_mem 20 -t $targetGenome -q $queryGenome -o $outDir

