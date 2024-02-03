#!/usr/bin/env bash
#SBATCH -J svim
#SBATCH --partition=long
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8

# sorted and indexed bam = $1
# ref assembly = $2

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate svim_env

fileshort=$(basename $1 .sorted.bam)

svim alignment $fileshort $1 $2
