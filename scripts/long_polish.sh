#!/usr/bin/env bash
#SBATCH -J long_polish
#SBATCH --partition=medium
#SBATCH --mem=180G
#SBATCH --cpus-per-task=8

Reads=$1
Assembly=$2

Short=$(basename $Assembly .fasta)

# Racon (x1)
/mnt/shared/scratch/jnprice/apps/raconnn/raconnn 1 $Reads $Assembly > "$Short"_racon.fasta


# Medaka (x1)
export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate medaka

mkdir $Short
medaka_consensus -i $Reads -d "$Short"_racon.fasta -o ./$Short -t 8 -m r941_min_high_g360
cp $Short/consensus.fasta "$Short"_medaka.fasta
rm -rf $Short

