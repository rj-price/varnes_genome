#!/usr/bin/env bash
#SBATCH -J orthofinder
#SBATCH --partition=medium
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

# proteome_dir = $1

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate orthofinder

ulimit -Sn
ulimit -Hn
ulimit -n 4096
ulimit -Sn

/mnt/shared/scratch/jnprice/apps/OrthoFinder/orthofinder.py -M msa -t 15 -a 30 -f $1 
