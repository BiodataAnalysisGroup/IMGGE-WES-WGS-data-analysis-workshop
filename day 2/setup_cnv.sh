#use environment.yml to create a conda environment named cnv
conda env create -f environment.yml

#activate environment
conda activate cnv

#check if cnv has been installed
cnvkit.py -h