# Day 1

# A. Quick overview of basic Unix commands (Navigating Files and Directories, Working With Files and Directories, Pipes and Filters, Loops, Finding Things)  

## Introduction

There are many different command-line interfaces. In this tutorial we work with WSL2 and Bash ("Bourne again shell"). Some of the many things one can do utilizing the command-line are:

- Navigate the computational environment
- View content
- Start and control the execution of programs
- Process files and folders
- Write and run code

## Navigating Files and Directories

Paths in linux can be "partial" (e.g., ../Documents) or "full" (e.g., /home/user/Documents).

Some useful commands for directory/file processing are:

- `ls`: list files
- `cd, pwd`: change directory, print working directory
- `cp, mv, rm`:  copy, move, remove files
- `du, df`:  look at disk usage
- `more, less`:  look at (text) files
- `find`: find files
- `grep`: find things with in files
- `mkdir, rmdir`:  make and remove directories
- `chmod`: change permissions
- `wget`: download content
- `tar`: unpack compressed content

## Working With Files and Directories

- `cat`: building datasets
- `sed`: reformatting data
- `awk`: powerful line by line processing of data
- `grep`: find patterns in data
- `head, tail`: look at tops and bottoms
- `sort, uniq`: arrange and simplify data
- `wc`: count elements in data
- `top`: see what processes are running
- `kill`: kills a given process number
- `python`: open Python interface
- `R`: open R interface

## Pipes and Filters, Loops, Finding Things

### Standard streams are input and output channels for a process

- `stdin` is the standard input stream (often the keyboard by default)
- `stdout` is the standard output stream (often the display by default)
- `stderr` is the standard error message stream (often the display by default)

### Input, Output and Direction

- `command1 > outfile` redirects the output of command1
to a file outfile.

- `command1 < infile` redirects the input for command1
to be from a file infile.

- `command1 < infile > outfile` combines the two.

- The "pipe" `command1 | command2`: redirects the input for
command2 to be from the output of command1.

### Examples

echo sends its argument to stdout, while cat sends stdin to stdout:

```bash
echo "Hello world!" # to screen
echo "Hello world!" > hw.txt # to file
cat < hw.txt # from file
echo "Hello world!" | cat # piped to cat, hence to screen
```
### If ... else loops

```bash
# set variable
file_of_interest=test1.txt
# check if file exists
if  [ -f $file_of_interest ]; then
    echo "${file_of_interest} exists."
else
    echo "${file_of_interest} does not exist."
fi
# set variable
dir_of_interest=dir1
# check if directory exists
if  [ -d $dir_of_interest ]; then
    echo "${dir_of_interest} exists."
else
    echo "${dir_of_interest} does not exist."
fi
```

### For loops

```bash
# iterate with for loop
for i in $(seq 1 3)
do
    # create new file
    touch file_${i}.test.txt
    # append to a single file
    echo $i >> list_of_numbers.txt
done
# print list of files
ls *.test.txt

# get number of lines
var1=$(cat list_of_numbers.txt | wc -l)
# iterate with for loop
for i in $(seq 1 ${var1})
do
    # print the contents of each line based on $i index
    echo $(sed "${i}q;d" list_of_numbers.txt)
done
```

# B. Structure of germline small variant calling workflow

## 1. List of software tools for installation (done automatically with conda and day1.yml)

- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trim galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) ([manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf))
- [BWA-MEM](https://anaconda.org/bioconda/bwa) ([manual](https://bio-bwa.sourceforge.net/bwa.shtml))
- [samtools](https://anaconda.org/bioconda/samtools) ([manual](https://www.htslib.org/doc/))
- [Picard tools](https://anaconda.org/bioconda/picard)
- [Freebayes](https://anaconda.org/bioconda/freebayes) ([manual](https://github.com/freebayes/freebayes))
- [bcftools](https://anaconda.org/bioconda/bcftools) ([manual](https://samtools.github.io/bcftools/bcftools.html))
- [SnpEff & SnpSift](http://pcingola.github.io/SnpEff/)
- GEMINI (load & inheritance pattern) through [Galaxy Europe](https://usegalaxy.eu/)

## 2. Background and Metadata (discussing the dataset to be used)

Using the Whole-Genome Sequencing (WGS) data from a family of three (father, mother and child) to run the tutorial on Variant calling workflow [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3243160.svg)](https://doi.org/10.5281/zenodo.3243160). The following tutorial is adjusted from the [Exome sequencing data analysis for diagnosing a genetic disease](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html) tutorial, making use of terminal commands rather than Galaxy tools, with the exception of [Gemini](https://doi.org/10.1371/journal.pcbi.1003153), to perform the analysis.

The files required for this tutorial can be downloaded directly from here:

- father ([R1](https://zenodo.org/record/3243160/files/father_R1.fq.gz?download=1)) ([R2](https://zenodo.org/record/3243160/files/father_R2.fq.gz?download=1))
- mother ([R1](https://zenodo.org/record/3243160/files/mother_R1.fq.gz?download=1)) ([R2](https://zenodo.org/record/3243160/files/mother_R2.fq.gz?download=1))
- child (proband) ([R1](https://zenodo.org/record/3243160/files/proband_R1.fq.gz?download=1)) ([R2](https://zenodo.org/record/3243160/files/proband_R2.fq.gz?download=1))
- [hg19_chr8.fa.gz](https://zenodo.org/record/3243160/files/hg19_chr8.fa.gz?download=1)
- [dbsnp_138.hg19.chr8.vcf.gz](https://zenodo.org/record/3243160/files/dbsnp_138.hg19.chr8.vcf.gz?download=1)
- [Pedigree.txt](https://zenodo.org/record/3243160/files/Pedigree.txt?download=1)

Pre-mapped BAM files are also available for all three samples (although mapping will be performed during the tutorial):

- [mapped_reads_father](https://zenodo.org/record/3243160/files/mapped_reads_father.bam?download=1)
- [mapped_reads_mother](https://zenodo.org/record/3243160/files/mapped_reads_mother.bam?download=1)
- [mapped_reads_child](https://zenodo.org/record/3243160/files/mapped_reads_proband.bam?download=1)

Download all files in the same folder, where the tutorial will take place. Decompress the files for reference genome (hg19_chr8.fa.gz) and dbSNP 138 (dbsnp_138.hg19.chr8.vcf.gz) using any appropriate tool (e.g., `gunzip`) beforehand to be ready for the tutorial (as it can be time-consuming).

## 3. Assessing Read Quality (discuss FASTQ format)

### FASTQC tool
```bash
### Introduction on Variant calling ###

# Different types of sequencing for different purposes. For example: 
# - Targeted Sequencing (e.g., COVID-19 and lineagespot?)
# - Whole-Exome Sequencing (WES)
# - Whole-Genome Sequencing (WGS)

# Goal of variant calling, i.e. identification of variants that are associated with the phenotype under study (e.g., genetic disorders, cancer)

### FastQC - Quality control ###

# Overview of FastQ format

# The need for quality control of raw data, following sequencing

# Running FastQC
fastqc --version
fastqc --help
fastqc *.fq.gz
```

## 4. Trimming and Filtering (discuss quality parameters and how to assess them)

### Trim Galore or Trimmomatic

```bash
### Trim galore / Trimmomatic - Trimming for low-quality, adapter and duplicated (e.g., optical or PCR duplicates, PCR bias) sequences ###
# Trim_galore
trim_galore --version
trim_galore --help
mkdir trim_galore_data
# set the option --path_to_cutadapt to the Cutadapt executable if Cutadapt is not in the PATH (default)
trim_galore --quality 25 -o trim_galore_data --gzip --fastqc --paired father_R1.fq.gz father_R2.fq.gz
trim_galore --quality 25 -o trim_galore_data --gzip --fastqc --paired mother_R1.fq.gz mother_R2.fq.gz
trim_galore --quality 25 -o trim_galore_data --gzip --fastqc --paired proband_R1.fq.gz proband_R2.fq.gz
# OR in a single command
trim_galore --quality 25 -o trim_galore_data --gzip --fastqc --paired \
father_R1.fq.gz father_R2.fq.gz \
mother_R1.fq.gz mother_R2.fq.gz \
proband_R1.fq.gz proband_R2.fq.gz

# Trimmomatic
mkdir trimmomatic_data
# paired-end mode
#
trimmomatic PE -threads 4 -phred33 father_R1.fq.gz father_R2.fq.gz \
trimmomatic_data/father_R1_paired.fq.gz trimmomatic_data/father_R1_unpaired.fq.gz \
trimmomatic_data/father_R2_paired.fq.gz trimmomatic_data/father_R2_unpaired.fq.gz \
ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
SLIDINGWINDOW:4:15 LEADING:5 TRAILING:3 MINLEN:36
#
trimmomatic PE -threads 4 -phred33 mother_R1.fq.gz mother_R2.fq.gz \
trimmomatic_data/mother_R1_paired.fq.gz trimmomatic_data/mother_R1_unpaired.fq.gz \
trimmomatic_data/mother_R2_paired.fq.gz trimmomatic_data/mother_R2_unpaired.fq.gz \
ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
SLIDINGWINDOW:4:15 LEADING:5 TRAILING:3 MINLEN:36
#
trimmomatic PE -threads 4 -phred33 proband_R1.fq.gz proband_R2.fq.gz \
trimmomatic_data/proband_R1_paired.fq.gz trimmomatic_data/proband_R1_unpaired.fq.gz \
trimmomatic_data/proband_R2_paired.fq.gz trimmomatic_data/proband_R2_unpaired.fq.gz \
ILLUMINACLIP:../Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
SLIDINGWINDOW:4:15 LEADING:5 TRAILING:3 MINLEN:36
```

## 5. Mapping and post-processing (discuss BAM/SAM format, mapping statistics and visualization of mapped reads)

### BWA-MEM, samtools, custom R script
```bash
### Mapping reads to reference genome - BWA-MEM ###

# Create index for reference genome
bwa index hg19_chr8.fa

# Map pre-processed reads with BWA-MEM
# Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments. (-a)
# Mark shorter split hits as secondary (for Picard compatibility). (-M)
bwa mem -t 4 -a -M hg19_chr8.fa trimmomatic_data/father_R1_paired.fq.gz trimmomatic_data/father_R2_paired.fq.gz > father.sam
bwa mem -t 4 -a -M hg19_chr8.fa trimmomatic_data/mother_R1_paired.fq.gz trimmomatic_data/mother_R2_paired.fq.gz > mother.sam
bwa mem -t 4 -a -M hg19_chr8.fa trimmomatic_data/proband_R1_paired.fq.gz trimmomatic_data/proband_R2_paired.fq.gz > proband.sam
```

```bash
### Pre-processing of mapped reads before variant calling ###

# Explain SAM & BAM format

# Convert SAM-to-BAM file
samtools view --threads 4 -b father.sam > father.bam
samtools view --threads 4 -b mother.sam > mother.bam
samtools view --threads 4 -b proband.sam > proband.bam

# Sort reads based on genomic coordinates
samtools sort --threads 4 -o father_sorted.bam father.bam
samtools sort --threads 4 -o mother_sorted.bam mother.bam
samtools sort --threads 4 -o proband_sorted.bam proband.bam

# Generate mapping statistics on BAM files
# 1. FASTQC on sorted BAM files
fastqc *_sorted.bam

# 2. samtools flagstat
samtools flagstat father_sorted.bam > father_sorted.flagstat
samtools flagstat mother_sorted.bam > mother_sorted.flagstat
samtools flagstat proband_sorted.bam > proband_sorted.flagstat
# OR (in case of large number of samples, use loops)
for i in $(ls *_sorted.bam)
do
    # extract sample name (remove suffix)
    sample=${i%_sorted*}
    # run command
    samtools flagstat $i > ${sample}.flagstat
done

# 3. samtools stats & plot-bamstat
for i in $(ls *_sorted.bam)
do
    # extract sample name (remove suffix)
    sample=${i%_sorted*} # keep the *_sorted* part of the filename
    # run samtools stat
    samtools stats $i > ${sample}.bc
    # extract summary statistics
    cat ${sample}.bc | grep ^SN | cut -f 2- > ${sample}.bc.summary
    # produce plot-bamstat plot
    plot-bamstats -p ./bamstat_plots/${sample} ${sample}.bc
done

# 4. custom script (e.g., in R) OR ready-to-use software solution (e.g., BAMStats (https://bamstats.sourceforge.net/))
Rscript summary_plot_generation.R

# Filter the paired-end reads of all samples to retain only those read pairs, 
# for which both the forward and the reverse read have been mapped to the reference successfully
# https://broadinstitute.github.io/picard/explain-flags.html

# exclude flags with -F
# require flags with -f
# read mapped in proper pair (0x2) - 2
# not primary alignment (0x100) - 256
samtools view --threads 4 -b -f 2 -o father_filtered.bam father_sorted.bam 
samtools view --threads 4 -b -f 2 -o mother_filtered.bam mother_sorted.bam
samtools view --threads 4 -b -f 2 -o proband_filtered.bam proband_sorted.bam
```

## 6. Variant calling (discuss VCF file format and annotation, Filtering / extraction of variants from defined genomic regions, impact of QC of the raw and mapped reads to the variant calling QC)

### Picard tools, Freebayes, bcftools

```bash
# Set read groups information (Picard)
picard AddOrReplaceReadGroups \
    I=father_filtered.bam \
    O=father_aor.bam \
    RGID=1 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=father
#
picard AddOrReplaceReadGroups \
    I=mother_filtered.bam \
    O=mother_aor.bam \
    RGID=2 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=mother
#
picard AddOrReplaceReadGroups \
    I=proband_filtered.bam \
    O=proband_aor.bam \
    RGID=3 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=proband

# MarkDuplicates (Picard)
picard MarkDuplicates \
      I=father_aor.bam \
      O=father_md.bam \
      M=father_md_metrics.txt
#
picard MarkDuplicates \
      I=mother_aor.bam \
      O=mother_md.bam \
      M=mother_md_metrics.txt
#
picard MarkDuplicates \
      I=proband_aor.bam \
      O=proband_md.bam \
      M=proband_md_metrics.txt

# Generate index (.fai) for the reference fasta
samtools faidx hg19_chr8.fa

# Generate indexes (.bai) for bam files
samtools index *_md.bam

### Variant calling and filtering - FreeBayes/GATK HaplotypeCaller/bcftools ###
# default settings (for low to high depth sequencing in haploid and diploid samples)
freebayes -f hg19_chr8.fa -C 2 -F 0.2 father_md.bam > father.vcf
freebayes -f hg19_chr8.fa -C 2 -F 0.2 mother_md.bam > mother.vcf
freebayes -f hg19_chr8.fa -C 2 -F 0.2 proband_md.bam > proband.vcf

# Variant filtering (bcftools) and index generation (tabix)
bcftools filter -s QUAL100 -e '%QUAL<100' -Oz father.vcf | bcftools filter -s DP5 -e 'DP<5' -m+ -Oz -o father_filtered.vcf.gz
tabix -p vcf father_filtered.vcf.gz
bcftools filter -s QUAL100 -e '%QUAL<100' -Oz mother.vcf | bcftools filter -s DP5 -e 'DP<5' -m+ -Oz -o mother_filtered.vcf.gz
tabix -p vcf mother_filtered.vcf.gz
bcftools filter -s QUAL100 -e '%QUAL<100' -Oz proband.vcf | bcftools filter -s DP5 -e 'DP<5' -m+ -Oz -o proband_filtered.vcf.gz
tabix -p vcf proband_filtered.vcf.gz

# Create list of VCF files
ls *_filtered.vcf.gz > tmp_vcf_file_list.txt

# Merge VCF 
bcftools merge --threads 4 --filter-logic + --missing-to-ref -Oz --file-list tmp_vcf_file_list.txt | bcftools view --threads 4 -f PASS -Ov -o merged.vcf

# Post-processing
bcftools norm --threads 4 -f hg19_chr8.fa --multiallelics -both -Ov -o merged_norm.vcf merged.vcf
```

## 7. Variant annotation

### SnpEff, SnpSift, Gemini (in Galaxy Europe)

```bash
# Variant annotation - SnpEff & SnpSift
# check available databases for Homo sapiens
snpEff databases | grep -i "Homo_sapiens"
# Build UCSC hg19 database
snpEff download -v hg19
# Annotate with snpEff hg19 genome and dbSNP138 - hg19 - chr8
snpEff ann -v -c snpEff/snpEff.config -noStats hg19 merged_norm.vcf | \
SnpSift Annotate -v -id dbsnp_138.hg19.chr8.vcf > merged_norm_anno.vcf
# Analysis of annotated variants with the 'Autosomal recessive' model of Gemini software
```
Run [Galaxy Europe tools](https://usegalaxy.eu/) `GEMINI load` and `GEMINI inheritance pattern` with the [tutorial-suggested parameters](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html#generating-a-gemini-database-of-variants-for-further-annotation-and-efficient-variant-queries)
