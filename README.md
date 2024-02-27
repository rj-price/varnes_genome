# Varnes Genome Analysis

Download ONT data from Novogene and check MD5. \
**All good!**

## Environmental variables
```
workDir=/mnt/shared/scratch/jnprice/varnes_genome
scriptsDir=/mnt/shared/scratch/jnprice/varnes_genome/scripts
```

## QC reads
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

Call SVs with [SVIM](https://github.com/eldariont/svim/wiki)
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

### QC assemblies
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

### Long Read Polishing
Polish scaffolded assembly with high quality reads.
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
    -query $workDir/longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta \
    -db mito \
    -max_target_seqs 3 -outfmt 6 -evalue 1e-3 -num_threads 2 \
    > mito_results.outfmt6

blastn \
    -query $workDir/longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta \
    -db plas1 \
    -max_target_seqs 3 -outfmt 6 -evalue 1e-3 -num_threads 2 \
    > plas1_results.outfmt6

blastn \
    -query $workDir/longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta \
    -db plas2 \
    -max_target_seqs 3 -outfmt 6 -evalue 1e-3 -num_threads 2 \
    > plas2_results.outfmt6

blastn \
    -query $workDir/longstitch/Rsp4_filt_necat_x80_purged_longstitch.fasta \
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

**Mitochondrial genome = ntLink_9, bctg00000119_1, bctg00000134_1, bctg00000154_1** \
**Plastid genome = bctg00000110_1**

### Dgenies
Align scaffolded assembly to MJ using Dgenies

**Chromosomal contigs** \
HiC_scaffold_1_RagTag   =   ntLink_0, ntLink_6 \
HiC_scaffold_2_RagTag   =   ntLink_1 \
HiC_scaffold_3_RagTag   =   ntLink_2 \
HiC_scaffold_4_RagTag   =   bctg00000004_1, ntLink_11, ntLink_3 \
HiC_scaffold_5_RagTag   =   ntLink_4, ntLink_10, ntLink_7,  \
HiC_scaffold_6_RagTag   =   bctg00000027_1, ntLink_5,  \
HiC_scaffold_7_RagTag   =   ntLink_8 

**Small repetitive contigs** \
bctg00000123_1  41,124bp \
bctg00000135_1  82,325bp \
bctg00000145_1  25,911bp \
bctg00000152_1  31,209bp

**Larger sequences found in other contigs** \
bctg00000095_1  457,930bp \
bctg00000111_1  164,386bp

### Polish with Illumina data

QC reads
```
for file in /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/*gz;
    do sbatch $scriptsDir/fastqc.sh $file
    done
```

Run Pilon
```
sbatch $scriptsDir/pilon_2lib-aa.sh \
    $workDir/long_polish/Rsp4_filt_necat_x80_purged_longstitch_medaka.fasta \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG55DSX3_L2_1.fq.gz \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG55DSX3_L2_2.fq.gz \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG72DSX3_L4_1.fq.gz \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG72DSX3_L4_2.fq.gz \
    ./ \
    3
```

## RepearModeler & RepeatMasker
Pull container & test
```
apptainer pull docker://dfam/tetools:1.88

apptainer exec --bind /mnt/shared:/mnt/shared -H /mnt/shared/home/jnprice $APPS/singularity_cache/tetools_1.88.sif RepeatModeler -help
```

Run RepeatModeller and RepeatMasker
```
sbatch $scriptsDir/repeatmasker.sh Rsp4_filt_necat_x80_purged_longstitch_medaka_pilon.fasta rubus
```

Calculate percentage repeats per contig
```
python $scriptsDir/repeats_per_contig.py Rsp4_filt_necat_x80_purged_longstitch_medaka_pilon.fasta.masked > repeats_per_contig.txt
```

bctg00000004_1_pilon: 55.34% \
bctg00000027_1_pilon: 20.85% \
bctg00000095_1_pilon: 58.57% \
bctg00000110_1_pilon: 19.77% \
bctg00000111_1_pilon: 65.60% \
bctg00000119_1_pilon: 6.86% \
bctg00000123_1_pilon: 13.91% \
bctg00000134_1_pilon: 6.02% \
bctg00000135_1_pilon: 38.99% \
bctg00000145_1_pilon: 30.64% \
bctg00000152_1_pilon: 11.96% \
bctg00000154_1_pilon: 27.62% \
ntLink_0_pilon: 14.31% \
ntLink_10_pilon: 63.04% \
ntLink_11_pilon: 28.83% \
ntLink_1_pilon: 46.56% \
ntLink_2_pilon: 46.70% \
ntLink_3_pilon: 16.71% \
ntLink_4_pilon: 17.39% \
ntLink_5_pilon: 46.38% \
ntLink_6_pilon: 54.63% \
ntLink_7_pilon: 37.93% \
ntLink_8_pilon: 50.06% \
ntLink_9_pilon: 8.93%


### Remove repetitive and organellar contigs from final genome file
Convert masked genome to singleline fasta
```
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Rsp4_filt_necat_x80_purged_longstitch_medaka_pilon.fasta.masked > Rsp4_filt_necat_x80_purged_longstitch_medaka_pilon.fasta.masked.sl
```

Extract genome contigs based on dgenies alignment, excluding organellar and repetitive contigs
```
cat contig.list | while read line; do 
    grep -w -A1 $line Rsp4_filt_necat_x80_purged_longstitch_medaka_pilon.fasta.masked.sl >> Rsp4_final.fasta
done
```

## Comparative Plots

### Dotplot
```
sbatch $scriptsDir/scripts/minimap.sh \
    $workDir/svim/RiMJ_ragtag_HiC.fasta \
    $workDir/final/Rsp4_final.fasta \
    $workDir/comparison_plots/dotplot

sbatch ~/scripts/shell_scripts/genome_comparison/dotplot/create_dotplot.sh \
    $workDir/comparison_plots/dotplot/minimap.paf \
    $workDir/comparison_plots/dotplot
```

### Circos Plot
```
sbatch $scriptsDir/satsuma_synteny.sh \
    $workDir/svim/RiMJ_ragtag_HiC.fasta \
    $workDir/final/Rsp4_final.fasta \
    $workDir/comparison_plots/circos/

# Generate chromosome length files
faidx $workDir/svim/RiMJ_ragtag_HiC.fasta -i chromsizes > RiMJ_names.csv
faidx $workDir/final/Rsp4_final.fasta -i chromsizes > Rsp4_names.csv

sbatch $scriptsDir/pyCircos.sh
```

## Annotation

### BRAKER3
Use Anitra RNA-seq data from AB/MJ paper and protein evidence from eudicot orthoDB combined with Rubus idaeus proteomes.

Index final assembly and align RNAseq reads using HISAT2.
```
hisat2-build \
    $workDir/final/Rsp4_final.fasta \
    $workDir/annotation/genome_index/Rsp4

for file in $workDir/annotation/rna_reads/*1.fq.gz; do
    short=$(realpath $file | sed 's/.R1_val_1.fq.gz//g')
    sbatch $scriptsDir/hisat2.sh $workDir/annotation/genome_index/Rsp4 $file $short.R2_val_2.fq.gz
done
```
Approx. 95% alignment rate


Run BRAKER3
```
sbatch $scriptsDir/BRAKER3.sh $workDir/final/Rsp4_final.fasta
```

QC proteome
```
sbatch $scriptsDir/busco_eudicots_proteome.sh Rsp4_prot.fasta


***** Results: *****

C:97.3%[S:78.6%,D:18.7%],F:0.3%,M:2.4%,n:2326      
2263    Complete BUSCOs (C)                        
1828    Complete and single-copy BUSCOs (S)        
435     Complete and duplicated BUSCOs (D)         
6       Fragmented BUSCOs (F)                      
57      Missing BUSCOs (M)                         
2326    Total BUSCO groups searched 

```

---
**TO DO:**
- Comparative plots (synteny)
- Put together report of results
---

<br>

# TOOL TESTING

## Quast Test
```
sbatch ../scripts/quast.sh \
    /mnt/shared/scratch/jnprice/varnes_genome/svim/RiMJ_ragtag_HiC.fasta \
    RiMJ_ragtag_HiC.gff3 \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG55DSX3_L2_1.fq.gz \
    /mnt/shared/projects/niab/dsargent/Rubus/Varnes_genome_sequence/X204SC22053791-Z01-F001/raw_data/R_Va/R_Va_EDSW220016932-1a_HKG55DSX3_L2_2.fq.gz \
    /mnt/shared/scratch/jnprice/varnes_genome/Rsp4_highqual.fastq.gz \
    . \
    /mnt/shared/scratch/jnprice/varnes_genome/final/Rsp4_final.fasta
```

## Synteny plots
```
# Generate chromosome length files
faidx ../svim/RiMJ_ragtag_HiC.fasta -i chromsizes > MJ_ls.len
faidx ../long_polish/highqual/Rsp4_highqual_necat_purged_racon.fasta -i chromsizes > Rsp4_ls.len

# Extract ntLink_4_pilon and HiC_scaffold_5_RagTag

# Reverse comp ntLink_4_pilon
seqkit seq -r -p -t DNA ntLink_4_pilon.fasta -o ntLink_4_pilon_rev.fasta

# Align
minimap2 -x asm20 -t 2 HiC_scaffold_5_RagTag.fasta ntLink_4_pilon_rev.fasta > minimap_asm20.paf

# Convert PAFs to link files
cat minimap_asm20.paf | cut -f1,3,4,6,8,9 > Rsp4_vs_MJ_asm20.link

# Run NGenomeSyn locally
```

