# use cnv_environment.yml to create a conda environment named cnv
conda env create -f cnv_environment.yml

# use plink_environment.yml to create a conda environment named plink_env
conda env create -f plink_environment.yml

# activate environment
conda activate cnv

# check if cnv has been installed
cnvkit.py -h

# activate environment
conda activate plink_env

# check if plink and r-qqman have been installed
plink -h
R --version
Rscript -e "library(qqman)"
