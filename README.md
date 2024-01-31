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
Identify and remove haplotigs and overlapping heterozygous contigs.
```
for file in $workDir/purge_dups/*fasta;
    do sbatch $scriptsDir/purge_dups.sh $file \
        $workDir/Rsp4_highqual.fastq.gz
    done
```

### Long Read Polishing
Polish purged assembly with filtered reads and high quality reads and compare.

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