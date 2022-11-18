# Serbia-WES-WGS-data-analysis
Training material for the WES/WGS data analysis course

# Installation notes for WSL2 and Docker Desktop

To install the computational environment required for the execution of the present tutorial, follow the instructions of the [cwltool installation tutorial](https://github.com/common-workflow-language/cwltool), specifically the section on **MS Windows users**. There you will find detailed instructions on how to activate and install:

- Windows Subsystem for Linux 2 (WSL2)
- Docker Desktop
- Ubuntu or Debian from the Microsoft Store

# Installation guidelines for software dependencies
Following the installation of WSL2, open a WSL instance and execute within the `/installation_scripts` directory the following:

```bash
bash wes_wgs_tutorial_software_installation.sh
```

to automatically set up and check Miniconda and the conda environments for: 

1. Day 1 
2. CNVkit (Day 2)
3. PLINK 1.9 (Day 2)

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

## Day 3

