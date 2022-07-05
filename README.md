# Abundance_determinants_snoRNA
__Author__ : Etienne Fafard-Couture

__Email__ :  _<etienne.fafard-couture@usherbrooke.ca>_

## Description

Snakemake-based workflow to predict the abundance status of human snoRNAs and identify their main abundance determinants. 
This pipeline also predicts by default the abundance status of mouse snoRNAs and can also used to predict the abundance 
status of several vertebrate species (see below).

## Requirements
* Conda (Tested with version=4.11.0)
* Mamba (Tested with version=0.15.3)
* Snakemake (Tested with version=7.1.0)

1 - Conda (Miniconda3) needs to be installed (https://docs.conda.io/en/latest/miniconda.html).
For Linux users:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Answer `yes` to `Do you wish the installer to initialize Miniconda3?`

2 - Mamba needs to be installed via conda (mamba greatly speeds up environment creation,
    but conda can still be used instead of mamba):
```bash
conda install -n base -c conda-forge mamba
```
3 - Snakemake needs to be installed via mamba:
```bash
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
To activate the 'snakemake' environment that was just created:
```bash
conda activate snakemake
```

## Running the pipeline on Slurm-based cluster (to predict the abundance status of human and mouse snoRNAs)
#### Firstly, download environments need to be created:
```bash
snakemake all_downloads --conda-create-envs-only --use-conda --conda-frontend mamba --cores 1
```
#### Secondly, datasets need to be downloaded:
```bash
snakemake all_downloads --use-conda --conda-frontend mamba --cores 1
```
#### Thirdly, environments required by tools and calculations need to be created:
```bash
snakemake --conda-create-envs-only --use-conda --conda-frontend mamba --cores 1
```
#### Fourthly, running the pipeline using computation nodes and the previously created envs is done as follows:
```bash
snakemake -j 999 --use-conda --immediate-submit --notemp --cluster-config cluster.json --cluster 'python3 slurmSubmit.py {dependencies}'
```
