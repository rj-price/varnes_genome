# Varnes Genome Analysis

Download ONT data from Novogene and check md5s. \
**All good!**

## Environmental variables
```
workDir=/mnt/shared/scratch/jnprice/varnes_genome
scriptsDir=/mnt/shared/scratch/jnprice/varnes_genome/scripts
```

## QC
Uncompress ```fastq_pass.tgz``` file. Trim adapters with Porechop and gzip trimmed reads.
```
sbatch $scriptsDir/porechop.sh \
    /mnt/shared/projects/niab/jnprice/rubus/Varnes/20240111_nanopore/ONT/Rsp4/TJPROJ6/TGS/haiwai/haiwai/HW_ONT_qc/X204SC23116571-Z01-F001/data_release/X204SC23116571-Z01-F001/raw_data/Rsp4/20231221_1538_4H_PAS75777_f65d910b/fastq_pass \
    $workDir/Rsp4_trim.fastq

gzip $workDir/Rsp4_trim.fastq
```

Check trimmed reads with NanoPlot
```
sbatch $scriptsDir/nanoplot_fastq.sh \
    $workDir/Rsp4_trim.fastq.gz \
    $workDir/nanoplot/trimmed/
```

### Filtering
Generate two read sets based on quality scores using Filtlong and QC using NanoPlot.

Q12+ and 10kb+
```
sbatch $scriptsDir/filtlong.sh $workDir/Rsp4_trim.fastq.gz

sbatch $scriptsDir/nanoplot_fastq.sh \
    $workDir/Rsp4_filt.fastq.gz \
    $workDir/nanoplot/filt/
```

Q17+ and 1kb+
```
sbatch $scriptsDir/filtlong_highqual.sh $workDir/Rsp4_trim.fastq.gz

sbatch $scriptsDir/nanoplot_fastq.sh \
    $workDir/Rsp4_highqual.fastq.gz \
    $workDir/nanoplot/highqual/
```

## Assembly

### NECAT
Perform assemblies on different read sets (Q12+ and Q17+), with different coverage (30x and 80x).

Generate config files
```
bash $scriptsDir/necat_config.sh $workDir/Rsp4_filt.fastq.gz
bash $scriptsDir/necat_config.sh $workDir/Rsp4_highqual.fastq.gz
```

Edit config files
||Q12+ 30x |Q12+ 80x|Q17+ 30x |Q17+ 80x |
|--|--|--|--|--|
|PROJECT=|Rsp4_filt|Rsp4_filt_x80 |Rsp4_highqual | Rsp4_highqual_x80 | |
|ONT_READ_LIST=|read_list.txt|read_list.txt |read_list.txt |read_list.txt | |
|GENOME_SIZE=|270000000|270000000 |270000000 |270000000 | |
|THREADS=|16|16 |16 |16 | |
|MIN_READ_LENGTH=|10000|10000 |10000 |10000 | |
|PREP_OUTPUT_COVERAGE=|40|80 |40 |80 | |
|CNS_OUTPUT_COVERAGE=|30|80 | 30|80 | |

<br>

Run assemblies
```
sbatch $scriptsDir/necat.sh config.txt
```

QC readsets used for each assembly
```
for folder in $workDir/assemblies/*/R*/R*/;
    do short=$(basename $folder)
    mkdir -p $workDir/nanoplot/assemblies/$short
    sbatch $scriptsDir/nanoplot_fasta.sh \
        $folder/trimReads.fasta.gz \
        $workDir/nanoplot/assemblies/$short
    done
```

### QC
```
for file in *.fasta; 
    do sbatch $scriptsDir/busco_eudicots.sh $file
    done

for file in *.fasta;
    do python $scriptsDir/assembly_stats.py $file
    done \
    && python $scriptsDir/assembly_stats_concat.py *stats

cp BUSCO_*/short* . && for file in *txt;
    do python $scriptsDir/busco_extract.py $file
    done \
    && python $scriptsDir/busco_concat.py *busco
```

### Purge Heterozygous Contigs
Identify and remove haplotype-specific assembly segments incorrectly labeled as primary contigs, as well as heterozygous contig overlaps.
```
for file in $workDir/purge_dups/*fasta;
    do sbatch $scriptsDir/purge_dups.sh $file \
        $workDir/Rsp4_highqual.fastq.gz
    done
```

### Long Read Polishing
Polish purged assemblies with filtered reads and high quality reads, and compare.
```
for file in $workDir/long_polish/*fasta;
    do sbatch $scriptsDir/long_polish.sh \
        $workDir/Rsp4_filt.fastq.gz $file
    done

for file in $workDir/long_polish/*fasta;
    do sbatch $scriptsDir/long_polish.sh \
        $workDir/Rsp4_highqual.fastq.gz $file
    done
```

### Polish with Illumina data

#### QC reads
```
for file in /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/*gz;
    do sbatch $scriptsDir/fastqc.sh $file
    done
```

---
# TOOL TESTING

#### Pilon test
```
sbatch ../scripts/pilon_2lib-aa.sh \
    ~/scratch/varnes_genome/long_polish/highqual/Rsp4_highqual_necat_purged_racon.fasta \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG55DSX3_L2_1.fq.gz \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG55DSX3_L2_2.fq.gz \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG72DSX3_L4_1.fq.gz \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG72DSX3_L4_2.fq.gz \
    ~/scratch/varnes_genome/pilon_test/ \
    3
```

## Quast Test
```
sbatch ../scripts/quast.sh \
    $workDir/svim/RiMJ_ragtag_HiC.fasta \
    RiMJ_ragtag_HiC.gff3 \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG55DSX3_L2_1.fq.gz \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG55DSX3_L2_2.fq.gz \
    $workDir/Rsp4_highqual.fastq.gz \
    . \
    /mnt/shared/scratch/jnprice/varnes_genome/long_polish/highqual/Rsp4_highqual_necat_purged_racon.fasta
```

## Structural Variants
Align reads using LRA, sort and index bam, and calculate 1Mb coverage
```
lra index -ONT RiMJ_ragtag_HiC.fasta

sbatch $scriptsDir/long_read_aligner.sh $workDir/Rsp4_highqual.fastq.gz $workDir/svim/RiMJ_ragtag_HiC.fasta
```

Call SVs with SVIM (https://github.com/eldariont/svim/wiki)
```
sbatch $scriptsDir/svim.sh $workDir/svim/Rsp4_highqual.sorted.bam $workDir/svim/RiMJ_ragtag_HiC.fasta
```

### Dotplot Test **WORKED**
```
sbatch minimap.sh ../svim/RiMJ_ragtag_HiC.fasta ../long_polish/highqual/Rsp4_highqual_necat_purged_racon.fasta .

sbatch ~/scripts/shell_scripts/genome_comparison/dotplot/create_dotplot.sh \
    /mnt/shared/scratch/jnprice/varnes_genome/dotplot_test/minimap.paf \
    /mnt/shared/scratch/jnprice/varnes_genome/dotplot_test
```

Sat  3 Feb 10:22:32 GMT 2024
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          15817832     himem long_pol  jnprice  R 1-02:39:17      1 n17-28-1536-apollo
          15817833     himem long_pol  jnprice  R 1-02:39:17      1 n17-28-1536-apollo
          15817835      long    quast  jnprice  R 1-02:38:47      1 n19-32-192-ultron
          15817834      long    pilon  jnprice  R 1-02:42:02      1 n19-32-192-taserface
          15826782      long     svim  jnprice  R       6:24      1 n19-32-192-mantis
          15826759    medium purge_du  jnprice  R      15:09      1 n19-32-192-ego