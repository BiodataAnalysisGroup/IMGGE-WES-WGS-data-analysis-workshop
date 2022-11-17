# Serbia-WES-WGS-data-analysis
Training material for the WES/WGS data analysis course

# Installation notes for WSL2 and miniconda

To install the computational environment required for the execution of the present tutorial, follow the instructions of the [cwltool installation tutorial](https://github.com/common-workflow-language/cwltool), specifically the section on **MS Windows users**. There you will find detailed instructions on how to activate and install:

- Windows Subsystem for Linux 2 (WSL2)
- Docker Desktop
- Ubuntu or Debian from the Microsoft Store

## Day 1 

### Conda environment installation

Create a conda environment, containing all software dependencies for the 1st day of the tutorial, by running the following command:

```bash
conda env create -f day1.yml
```
The **day1.yml** is available in the `/day 1` folder.

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

## Day 3

