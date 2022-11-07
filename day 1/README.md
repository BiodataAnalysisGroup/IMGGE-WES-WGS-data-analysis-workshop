# A. Quick overview of basic Unix commands (Navigating Files and Directories, Working With Files and Directories, Pipes and Filters, Loops, Finding Things)  

# B. Structure of germline small variant calling

## 1. Background and Metadata (discussing the dataset to be used)

Using the Whole-Genome Sequencing (WGS) data from a family of three (father, mother and child) to run the tutorial on Variant calling workflow [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3243160.svg)](https://doi.org/10.5281/zenodo.3243160). The following tutorial is adjusted from the [Exome sequencing data analysis for diagnosing a genetic disease](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html) tutorial, making use of command-line commands, rather than Galaxy tools with the exception of [Gemini](https://doi.org/10.1371/journal.pcbi.1003153), to perform the analysis.

## 2. Assessing Read Quality (discuss FASTQ format)

FastQC tool

## 3. Trimming and Filtering (discuss quality parameters and how to assess them)

Trim Galore or Trimmomatic

## 4. Mapping and post-processing (discuss BAM/SAM format, mapping statistics and visualization of mapped reads)

BWA-MEM, samtools, Picard tools

## 5. Variant calling (discuss VCF file format and annotation, Filtering / extraction of variants from defined genomic regions, impact of QC of the raw and mapped reads to the variant calling QC)

Freebayes, bcftools

## 6. Variant annotation

SnpEff, SnpSift, Gemini (in Galaxy Europe)
