#!/usr/bin/env bash
#SBATCH -J busco
#SBATCH --partition=medium
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

# assembly (in .fasta format) = $1

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate busco

fileshort=$(basename $1 | sed s/".fasta"//g)
busco -m genome -c 8 -i $1 -o BUSCO_$fileshort -l /mnt/shared/home/jnprice/busco_downloads/lineages/eudicots_odb10

