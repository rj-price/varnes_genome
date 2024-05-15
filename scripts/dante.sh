#!/bin/bash
#SBATCH -J dante
#SBATCH --partition=medium
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8

# INPUTS
Assembly=$1

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate dante

OutDir=$(basename $1 .fasta)

~/scratch/apps/dante/dante.py \
    --query $Assembly \
    --protein_database /mnt/shared/scratch/jnprice/apps/dante/tool-data/protein_domains/Viridiplantae_v3.0_pdb \
    --classification /mnt/shared/scratch/jnprice/apps/dante/tool-data/protein_domains/Viridiplantae_v3.0_class \
    --output_dir $OutDir
