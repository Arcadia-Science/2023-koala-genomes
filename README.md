# Profiling expansions in [koala population seq data](https://github.com/awgg-lab/australasiangenomes/blob/main/species/Phascolarctos_cinereus.md) 🐨
The pipeline described in this repository was created to profile the polyglutamine and polyalanine repeat expansion lengths in the RUNX2, ZIC2, FOXL2, and ARX genes across a sequenced koala population. It is associated with the [Repeat expansions associated with human disease are present in diverse organisms](https://doi.org/10.57844/arcadia-e367-8b55) pub from Arcadia Science.

### Set-up
#### First create a conda environment using the ```envs/environment.yml``` file and activate your new environment
This repository uses snakemake to run the pipeline and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html). After installing conda, run the following commands to create and activate an environment for running this pipeline:
```
conda env create -f envs/environment.yml -n koala
conda activate koala
```
#### Choose what bam files you want to download from the koala pop seq AWS bucket
The ```inputs/samples.tsv``` file contains all of the necessary information for the samples that will go through this workflow. If you want to run this workflow on a different set of files, you will need to edit this table. Each bam file of interest should have a single row in the table. These two columns are required for running the pipeline:

```AWSFileName```: File name from this [table](https://koalagenomes.s3.ap-southeast-2.amazonaws.com/Koala_Metadata-19-10-2022.csv). <br>
```AWSFolderName```: Corresponding folder from this [table](https://koalagenomes.s3.ap-southeast-2.amazonaws.com/Koala_Metadata-19-10-2022.csv). <br>


#### Run Snakemake from your activated environment
```
snakemake -p -j 8 --use-conda
```

```-p``` tells Snakemake to print the shell commands that it runs for each rule, which can be useful for understanding what is happening/debugging.<br>
```-j``` tells Snakemake how many cores to use for executing jobs in parallel. Adjust to your system. <br>
```--use-conda``` tells Snakemake to utilize Conda environments to manage dependencies for individual rules.<br>
The ```--rerun-incomplete``` flag can be appended to the above command if the pipeline terminates prematurely for some reason and needs to be restarted.<br>

#### Expected final output
A table of the polyA expansion lengths for RUNX2, ZIC2, FOXL2, and ARX. PolyQ expansion is also counted for RUNX2.

|sample|  length|  gene_aminoacid|
| ----------- | ----------- |----------- |
|Kyogle_F_50433|  22|      RUNX2_a|
|Kyogle_F_50433|  18|      ZIC2_a|
|Kyogle_F_50433|  16|      FOXL2_a|
|Kyogle_F_50433|  19|      ARX_a|
