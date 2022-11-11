# A. Quick overview of basic Unix commands (Navigating Files and Directories, Working With Files and Directories, Pipes and Filters, Loops, Finding Things)  

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

Using the Whole-Genome Sequencing (WGS) data from a family of three (father, mother and child) to run the tutorial on Variant calling workflow [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3243160.svg)](https://doi.org/10.5281/zenodo.3243160). The following tutorial is adjusted from the [Exome sequencing data analysis for diagnosing a genetic disease](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html) tutorial, making use of command-line commands, rather than Galaxy tools with the exception of [Gemini](https://doi.org/10.1371/journal.pcbi.1003153), to perform the analysis.

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
