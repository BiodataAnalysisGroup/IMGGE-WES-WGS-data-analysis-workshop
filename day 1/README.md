# A. Quick overview of basic Unix commands (Navigating Files and Directories, Working With Files and Directories, Pipes and Filters, Loops, Finding Things)  

## Introduction

There are many different command-line interfaces. In this tutorial we work with Bash ("Bourne again shell"). Some of the many things one can do utilizing command-line are:

- Navigate computer
- View content
- Start and control the running of programs
- Process files
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

# B. Structure of germline small variant calling

## 1. List of software tools for installation

- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trim galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) ([manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf))
- [BWA-MEM](https://anaconda.org/bioconda/bwa) ([manual](https://bio-bwa.sourceforge.net/bwa.shtml))
- [samtools](https://anaconda.org/bioconda/samtools) ([manual](https://www.htslib.org/doc/))
- [Picard tools](https://anaconda.org/bioconda/picard)
- [Freebayes](https://anaconda.org/bioconda/freebayes) ([manual](https://github.com/freebayes/freebayes))
- [bcftools](https://anaconda.org/bioconda/bcftools) ([manual](https://samtools.github.io/bcftools/bcftools.html))
- [SnpEff & SnpSift](http://pcingola.github.io/SnpEff/) ([manual]())
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

## 3. Assessing Read Quality (discuss FASTQ format)

FastQC tool

## 4. Trimming and Filtering (discuss quality parameters and how to assess them)

Trim Galore or Trimmomatic

## 5. Mapping and post-processing (discuss BAM/SAM format, mapping statistics and visualization of mapped reads)

BWA-MEM, samtools, custom R script 

## 6. Variant calling (discuss VCF file format and annotation, Filtering / extraction of variants from defined genomic regions, impact of QC of the raw and mapped reads to the variant calling QC)

Picard tools, Freebayes, bcftools

## 7. Variant annotation

SnpEff, SnpSift, Gemini (in Galaxy Europe)
