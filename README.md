# WES/WGS data analysis
Training material for the WES/WGS data analysis course

# Installation notes for WSL2 and Docker Desktop

To install the computational environment required for the execution of the present tutorial, follow the instructions of the [cwltool installation tutorial](https://github.com/common-workflow-language/cwltool), specifically the section on **MS Windows users**. There you will find detailed instructions on how to activate and install:

- Windows Subsystem for Linux 2 (WSL2)
- Docker Desktop
- Ubuntu or Debian from the Microsoft Store

# Installation guidelines for software dependencies
Following the installation of WSL2, open a WSL instance within the `/installation_scripts` directory and execute the following:

```bash
bash wes_wgs_tutorial_software_installation.sh
```

to automatically set up and check Miniconda and the conda environments for: 

1. **Variant calling workflow** (Day 1)
2. **CNVkit** (Day 2)
3. **PLINK 1.9** (Day 2)

The `/installation_scripts` directory contains the required YAML files for this purpose.

## Day 1 

### Required files

The files required for day1 of the tutorial can be downloaded directly from here:

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

## Day 2

### Required files

Some of the files required for day2 of the tutorial can be downloaded directly from here:

- [refFlat.txt](https://github.com/BiodataAnalysisGroup/IMGGE-WES-WGS-data-analysis-workshop/blob/main/day%202/useful_files/refFlat.txt)

- [access regions .bed file](https://github.com/BiodataAnalysisGroup/IMGGE-WES-WGS-data-analysis-workshop/blob/main/day%202/useful_files/access-5kb-mappable.hg19_chr5_chr12_chr17.bed)

- [BED file for the captured reagions](https://zenodo.org/record/5697358/files/capture_targets_chr5_12_17.bed)

The original `fasta/fastq` files can be found here:

- [reference file](https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz)

Samples:

- [Normal_R1](https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz)

- [Normal_R2](https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz)

- [Tumor_R1](https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz)

- [Tumor_R2](https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz)

but the processed files, along with the corresponding VCF file, will be given to you during the workshop.

## Day 3

### Software dependencies

- A WSL2 distro (Ubuntu/Debian) 
- Docker Desktop (connected to WSL2)
- A CWL runner (e.g., cwltool)

### Overview

Reproducible research is an important part of good scientific practice. Establishing robust computational analysis workflow facilitates reproducible research. There has been a strong emphasis on establishing graphic workflow management systems, such as Galaxy and Knime. However, the landscape and complexity of Linux command line tools for workflow management is vast and complicated. We aim to take advantage of European/International leaders in the field to present current scientific workflow paradigms and a lead workshops to implement a basic workflow. Topics covered include workflows and containerization.

### Learning objectives

### Docker
- Understand what containerisation is, and why you might use it in bioinformatics.
- Be familiar with Docker; basic concepts and structure.
- Find and run containers built by other people.

### CWL
- Understand syntax and structure for CWL.
- Understand how to write CWL tool definitions for command line tools.
- Read and write CWL files written in YAML.
- Join CWL tools into a workflow.
- Use Docker with CWL to provide software dependencies and ensure reproducibility.

### Audience and requirements

This introductory tutorial is aimed towards bioinformaticians (graduate students and researchers), who are interested in becoming familiar to Docker based workflows in CWL.

### Prerequisites

- Experience in Shell; this includes basic commands (such as `ls`, `cp`, `mv`, `nano/vim`) and operations such as (`apt`, installing tools etc).

**Given that the workshop will offer hands-on training, participants will need to bring their own laptop with them.**
