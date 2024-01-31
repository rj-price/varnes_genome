# Varnes Genome Analysis

Download ONT data from Novogene and check md5s. **All good!**

## QC
Uncompress ```fastq_pass.tgz``` file. Trim adapters with Porechop and gzip trimmed reads.
```
sbatch /mnt/shared/scratch/jnprice/rubus/porechop.sh \
    /mnt/shared/projects/niab/jnprice/rubus/Varnes/20240111_nanopore/ONT/Rsp4/TJPROJ6/TGS/haiwai/haiwai/HW_ONT_qc/X204SC23116571-Z01-F001/data_release/X204SC23116571-Z01-F001/raw_data/Rsp4/20231221_1538_4H_PAS75777_f65d910b/fastq_pass \
    /mnt/shared/scratch/jnprice/rubus/Rsp4_trim.fastq

gzip Rsp4_trim.fastq
```

Check trimmed reads with NanoPlot
```
sbatch /home/jnprice/scripts/shell_scripts/read_qc/nanoplot_fastq.sh \
    /mnt/shared/scratch/jnprice/rubus/Rsp4_trim.fastq.gz \
    /mnt/shared/scratch/jnprice/rubus/nanoplot/trimmed/
```

Remove reads <10kb and Q<12 with Filtlong
```
sbatch filtlong.sh /mnt/shared/scratch/jnprice/rubus/Rsp4_trim.fastq.gz
```

Check filtered reads with NanoPlot
```
sbatch /home/jnprice/scripts/shell_scripts/read_qc/nanoplot_fastq.sh \
    /mnt/shared/scratch/jnprice/rubus/Rsp4_filt.fastq.gz \
    /mnt/shared/scratch/jnprice/rubus/nanoplot/filt/
```

## Assembly

### NECAT
Generate config file
```
bash ~/scripts/shell_scripts/genome_assembly/necat_config.sh /mnt/shared/scratch/jnprice/rubus/Rsp4_filt.fastq.gz
```

Edit config file
```
PROJECT=Rsp4_filt
ONT_READ_LIST=read_list.txt
GENOME_SIZE=270000000
THREADS=16
```

Run assembly
```
sbatch ~/scripts/shell_scripts/genome_assembly/necat.sh Rsp4_filt_config.txt
```

### QC
```
for file in *.fasta; 
    do sbatch /home/jnprice/scripts/shell_scripts/genome_qc/busco_eudicots.sh $file
    done

for file in *.fasta;
    do python /home/jnprice/scripts/shell_scripts/genome_qc/assembly_stats.py $file
    done && python /home/jnprice/scripts/shell_scripts/genome_qc/assembly_stats_concat.py *stats

cp BUSCO_*/short* . && for file in *txt;
    do python /home/jnprice/scripts/shell_scripts/genome_qc/busco_extract.py $file
    done && python /home/jnprice/scripts/shell_scripts/genome_qc/busco_concat.py *busco
```

### Purge Heterozygous Contigs

```
for file in /mnt/shared/scratch/jnprice/rubus/purge_dups/*fasta;
    do sbatch /home/jnprice/scripts/shell_scripts/genome_assembly/purge_dups.sh $file \
        /mnt/shared/scratch/jnprice/rubus/Rsp4_highqual.fastq.gz
    done
```


## Polishing

Select high quality reads (>1kb and >Q17) for polishing
```
sbatch filtlong_highqual.sh /mnt/shared/scratch/jnprice/rubus/Rsp4_trim.fastq.gz
```

Check high quality reads with NanoPlot
```
sbatch /home/jnprice/scripts/shell_scripts/read_qc/nanoplot_fastq.sh \
    /mnt/shared/scratch/jnprice/rubus/Rsp4_highqual.fastq.gz \
    /mnt/shared/scratch/jnprice/rubus/nanoplot/highqual/
```


### Long Read Polishing
Polish purged assembly with filtered reads and high quality reads and compare

```
for file in /mnt/shared/scratch/jnprice/rubus/long_polish/*fasta;
    do sbatch /home/jnprice/scripts/shell_scripts/genome_assembly/long_polish.sh \
        /mnt/shared/scratch/jnprice/rubus/Rsp4_filt.fastq.gz $file
    done

for file in /mnt/shared/scratch/jnprice/rubus/long_polish/*fasta;
    do sbatch /home/jnprice/scripts/shell_scripts/genome_assembly/long_polish.sh \
        /mnt/shared/scratch/jnprice/rubus/Rsp4_highqual.fastq.gz $file
    done
```