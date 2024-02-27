#!/usr/bin/env bash
#SBATCH -J minimap
#SBATCH --partition=medium
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4

#INPUT:
targetGenome=$1
queryGenome=$2
outDir=$3

minimap2 -x asm5 -t 8 $targetGenome $queryGenome > $outDir/minimap.paf
