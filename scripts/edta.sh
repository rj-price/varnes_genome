#!/bin/bash
#SBATCH -J edta
#SBATCH --partition=long
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8

# INPUTS
Assembly=$(basename $1)
CDS=$2

export MYCONDAPATH=/mnt/shared/scratch/jnprice/apps/conda
source ${MYCONDAPATH}/bin/activate EDTA2

#apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/EDTA.sif 

~/scratch/apps/EDTA/EDTA.pl \
    --genome $Assembly \
    --cds $CDS \
    --threads 15 \
    --sensitive 1 \
    --anno 1
