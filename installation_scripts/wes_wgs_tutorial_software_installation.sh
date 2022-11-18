# 1. Installation of WSL2 and Docker Desktop

# To install the computational environment required for the execution of the present tutorial, follow the instructions of the cwltool installation tutorial, specifically the section on MS Windows users. 
# There you will find detailed instructions on how to activate and install:

# https://github.com/common-workflow-language/cwltool

# 1.1) Windows Subsystem for Linux 2 (WSL2)
# 1.2) Docker Desktop
# 1.3) Ubuntu or Debian from the Microsoft Store

# 2. Miniconda installation in WSL2

# Run the following commands in your WSL distro:
sudo apt-get update
sudo apt-get install wget
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh

# Verify shell script download using the following command line:
sha256sum Miniconda3-py39_4.12.0-Linux-x86_64.sh
# 78f39f9bae971ec1ae7969f0516017f2413f17796670f7040725dd83fcff5689  Miniconda3-py39_4.12.0-Linux-x86_64.sh

# Install Miniconda
sh ./Miniconda3-py39_4.12.0-Linux-x86_64.sh

# 3. Day 1 conda environment installation and validation
conda env create -f day1.yml
conda activate day1

# Validate that tools are functional
fastqc --help
trim_galore --help
trimmomatic --help
bwa
samtools --help
picard -h
bcftools --help
freebayes --help
snpEff --help
SnpSift --help

# 4. Day 2 conda environment installation and validation

# 4.1) CNVkit conda environment
conda env create -f cnv_environment.yml
conda activate cnv

# Validate that tools are functional
cnvkit.py -h

# 4.2) PLINK 1.9 conda environment
conda env create -f plink_environment.yml
conda activate plink_env

# Validate that tools are functional
plink -h
R --version
Rscript -e "library(qqman)"
