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

## Analysis of Structural Variants
Align reads to Malling Jewel using LRA, sort and index bam, and calculate 1Mb coverage.
```
lra index -ONT RiMJ_ragtag_HiC.fasta

sbatch $scriptsDir/long_read_aligner.sh $workDir/Rsp4_highqual.fastq.gz $workDir/svim/RiMJ_ragtag_HiC.fasta
```

Call SVs with SVIM (https://github.com/eldariont/svim/wiki)
```
sbatch $scriptsDir/svim.sh $workDir/svim/Rsp4_highqual.sorted.bam $workDir/svim/RiMJ_ragtag_HiC.fasta
```

Filter variants
```
conda activate sam_bcf_tools_env

bcftools view variants.vcf -H -i 'QUAL>40 && SUPPORT>10 && FILTER="PASS"' | wc -l
# 42269

bcftools view variants.vcf -i 'QUAL>40 && SUPPORT>10 && FILTER="PASS"' > sv_filt.vcf
```

Insertion in ANS gene!!
```
HiC_scaffold_5_RagTag	7906987	svim.INS.46403	N	<INS>	44	PASS	SVTYPE=INS;END=7906987;SVLEN=4337;SUPPORT=36;STD_SPAN=10.63;STD_POS=5.44	GT:DP:AD	0/1:60:24,36
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

### Scaffolding using LongStitch
Pull & test container
```
apptainer pull docker://ebird013/longstitch:1.0.5

apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/longstitch_1.0.5.sif longstitch run
```

```
sbatch $scriptsDir/longstitch.sh \
    $workDir/purge_dups/Rsp4_filt_necat_x80_purged.fasta \
    $workDir/Rsp4_filt.fastq.gz \
    270m
```

### Dgenies
Align scaffolded assembly to MJ using Dgenies

HiC_scaffold_1_RagTag   =   ntLink_0, ntLink_6
HiC_scaffold_2_RagTag   =   ntLink_1
HiC_scaffold_3_RagTag   =   ntLink_2
HiC_scaffold_4_RagTag   =   bctg00000004_1, ntLink_11, ntLink_3
HiC_scaffold_5_RagTag   =   ntLink_4, ntLink_10, ntLink_7, 
HiC_scaffold_6_RagTag   =   bctg00000027_1, ntLink_5, 
HiC_scaffold_7_RagTag   =   ntLink_8

### Long Read Polishing
Polish purged assemblies with filtered reads and high quality reads, and compare.
```
sbatch $scriptsDir/long_polish.sh \
    $workDir/Rsp4_highqual.fastq.gz \
    $workDir/longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta
```

### Organelle Check
Download NCBI mitochondial and chloroplast databases
```
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.3.1.genomic.fna.gz
```
Uncompress and make BLAST dbs
```
makeblastdb -in mitochondrion.1.1.genomic.fna -dbtype nucl -out mito
makeblastdb -in plastid.1.1.genomic.fna -dbtype nucl -out plas1
makeblastdb -in plastid.2.1.genomic.fna -dbtype nucl -out plas2
makeblastdb -in plastid.3.1.genomic.fna -dbtype nucl -out plas3
```

BLASTn contigs to organelle dbs
```
blastn \
    -query ../longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta \
    -db mito \
    -max_target_seqs 3 -outfmt 6 -evalue 1e-3 -num_threads 2 \
    > mito_results.outfmt6

blastn \
    -query ../longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta \
    -db plas1 \
    -max_target_seqs 3 -outfmt 6 -evalue 1e-3 -num_threads 2 \
    > plas1_results.outfmt6

blastn \
    -query ../longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta \
    -db plas2 \
    -max_target_seqs 3 -outfmt 6 -evalue 1e-3 -num_threads 2 \
    > plas2_results.outfmt6

blastn \
    -query ../longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta \
    -db plas3 \
    -max_target_seqs 3 -outfmt 6 -evalue 1e-3 -num_threads 2 \
    > plas3_results.outfmt6
```

Extract number of hits and total hit length
```
awk -F'\t' '{ count[$1]++; sum[$1]+=$4 } END { for (value in count) print value, count[value], sum[value] }' mito_results.outfmt6
awk -F'\t' '{ count[$1]++; sum[$1]+=$4 } END { for (value in count) print value, count[value], sum[value] }' plas1_results.outfmt6
awk -F'\t' '{ count[$1]++; sum[$1]+=$4 } END { for (value in count) print value, count[value], sum[value] }' plas2_results.outfmt6
awk -F'\t' '{ count[$1]++; sum[$1]+=$4 } END { for (value in count) print value, count[value], sum[value] }' plas3_results.outfmt6
```

**Mitochondrial genome = ntLink_9**
**Plastid genome = bctg00000110_1**

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


### Dotplot Test **WORKED**
```
sbatch minimap.sh ../svim/RiMJ_ragtag_HiC.fasta ../longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta .

sbatch ~/scripts/shell_scripts/genome_comparison/dotplot/create_dotplot.sh \
    /mnt/shared/scratch/jnprice/varnes_genome/dotplot_test/minimap.paf \
    /mnt/shared/scratch/jnprice/varnes_genome/dotplot_test
```

### Synteny plots
```
# Generate chromosome length files
faidx ../svim/RiMJ_ragtag_HiC.fasta -i chromsizes > MJ_ls.len
faidx ../long_polish/highqual/Rsp4_highqual_necat_purged_racon.fasta -i chromsizes > Rsp4_ls.len

# Remove alignments under 10kb
awk -F'\t' '$11 >= 10000' minimap.paf > minimap_10kb.paf

# Convert PAFs to link files
cat minimap_10kb.paf | cut -f1,3,4,6,8,9 > MJ_vs_Rsp4.link

# Run NGenomeSyn locally
```

## RepearModeler & RepeatMasker
Pull container & test
```
apptainer pull docker://dfam/tetools:1.88

apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/tetools_1.88.sif RepeatModeler -help
```

```
#! -pa = threads

BuildDatabase -name arabidopsis TAIR10_chr_all.fas

RepeatModeler -database arabidopsis -pa 16 -LTRStruct > out.log

RepeatMasker -pa 16 -gff -lib consensi.fa.classified -dir MaskerOutput TAIR10_chr_all.fas

```
