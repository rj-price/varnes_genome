#!/bin/bash
#SBATCH -J repeatmasker
#SBATCH --partition=long
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

# INPUTS
Assembly=$(basename $1)
db_name=$2

apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/tetools_1.88.sif \
    BuildDatabase -name $db_name $Assembly

apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/tetools_1.88.sif \
    RepeatModeler -database $db_name -threads 16 -LTRStruct > out.log

apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/tetools_1.88.sif \
    RepeatMasker -pa 16 -gff -lib rubus-families.fa -xsmall -html -gff -dir MaskerOutput $Assembly
