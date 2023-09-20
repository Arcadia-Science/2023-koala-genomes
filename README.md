# Profiling expansions in koala population seq data

### Set-up
#### First create a conda environment using the ```envs/environment.yml``` file and activate your new environment
This repository uses snakemake to run the pipeline and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html). After installing conda, run the following commands to create and activate an environment for running this pipeline:
```
conda env create -f envs/environment.yml -n chelicerate
conda activate chelicerate
```
