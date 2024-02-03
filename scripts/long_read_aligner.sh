#!/usr/bin/env bash
#SBATCH -J lra
#SBATCH --partition=medium
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8

# reads = $1
# ref assembly = $2

fileshort=$(basename $1 .fastq.gz)

zcat $1 | lra align -ONT -t 8 $2 /dev/stdin -p s > $fileshort.sam

samtools view -@ 8 -bS $fileshort.sam -o $fileshort.bam
samtools sort -@ 8 $fileshort.bam -o $fileshort.sorted.bam
samtools index -@ 8 $fileshort.sorted.bam

samtools depth $fileshort.sorted.bam > $fileshort.depth
python ~/scripts/my_python/average_window_depth.py $fileshort.depth $fileshort.1Mb_depth 1000000