# Serbia-WES-WGS-data-analysis
Training material for the WES/WGS data analysis course

# Installation notes for WSL2 and miniconda

To install the computational environment required for the execution of the present tutorial, follow the instructions of the [cwltool installation tutorial](https://github.com/common-workflow-language/cwltool), specifically the section on **MS Windows users**. There you will find detailed instructions on how to activate and install:

- Windows Subsystem for Linux 2 (WSL2)
- Docker Desktop
- Ubuntu or Debian from the Microsoft Store

## Day 1 conda environment installation

Create a conda environment, containing all software dependencies for the 1st day of the tutorial, by running the following command:

```bash
conda env create -f day1.yml
```
The **day1.yml** is available in the `/day 1` folder.
