#!/bin/bash
#SBATCH -J BRAKER3
#SBATCH --partition=long
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jordan.price@niab.com

Genome=$(realpath $1)
Prefix=$(basename $1 _final.fasta)
mkdir -p /mnt/shared/scratch/jnprice/varnes_genome/annotation/BRAKER3/$Prefix

RNA=$(ls /mnt/shared/scratch/jnprice/varnes_genome/annotation/alignment/*bam | tr '\n' ',')

apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/braker3_latest.sif braker.pl \
    --Species=rubus --gff3 --threads=16 \
	--AUGUSTUS_CONFIG_PATH=/mnt/shared/scratch/jnprice/apps/Augustus/config \
	--genome=$Genome \
	--bam=$RNA \
	--prot_seq=/mnt/shared/scratch/jnprice/varnes_genome/annotation/prot_db/eudicot_orthoDB_plusRubus_renamed2.fasta \
    --workingdir=/mnt/shared/scratch/jnprice/varnes_genome/annotation/BRAKER3/$Prefix