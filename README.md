# Abundance_determinants_snoRNA
__Author__ : Etienne Fafard-Couture

__Email__ :  _<etienne.fafard-couture@usherbrooke.ca>_

## Description

Snakemake-based workflow to predict the abundance status of human snoRNAs and identify their main abundance determinants. 
This pipeline also predicts by default the abundance status of mouse snoRNAs and can also used to predict the abundance 
status of several vertebrate species (see below).

## Requirements
* Conda (Tested with version=4.12.0)
* Mamba (Tested with version=0.15.3)
* Snakemake (Tested with version=6.0.5)

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
mamba create -c conda-forge -c bioconda -n snakemake snakemake=6.0.5
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
#### Secondly, datasets need to be downloaded (this might take a while):
```bash
snakemake all_downloads --use-conda --cores 1
```
#### Thirdly, environments required by tools and calculations need to be created:
```bash
snakemake --conda-create-envs-only --use-conda --conda-frontend mamba --cores 1
```
#### Fourthly, running the pipeline using computation nodes and the previously created envs is done as follows (Human snoRNA prediction):
```bash
snakemake -j 999 --use-conda --immediate-submit --notemp --cluster-config cluster.json --cluster 'python3 slurmSubmit.py {dependencies}'
```
#### Fifthly, running the pipeline using computation nodes and the previously created envs is done as follows (Mouse snoRNA prediction):                                                                           
```bash
snakemake all_mouse -j 999 --use-conda --immediate-submit --notemp --cluster-config cluster.json --cluster 'python3 slurmSubmit.py {dependencies}'
``` 
#### Finally, generating figures for the human and mouse snoRNAs is done as follows on the cluster:
```bash
snakemake all_figures -j 999 --use-conda --immediate-submit --notemp --cluster-config cluster.json --cluster 'python3 slurmSubmit.py {dependencies}'
```
#### or locally (of note: mouse figures might differ slightly from those in the paper since RNAcentral data is updated frequently):
```bash
snakemake all_figures --use-conda --cores 1
```


## Running the pipeline on Slurm-based cluster (to predict the abundance status of vertebrate species snoRNAs)
#### Firstly, the human/mouse snoRNA pipeline described above needs to be run (to generate the model used to predict species snoRNA abundance status)
#### Secondly, datasets need to be downloaded (this might take a while):
```bash
snakemake species_downloads --use-conda --cores 1
```
#### Thirdly, running the species prediction pipeline using computation nodes and the previously created envs is done as follows:
```bash
snakemake species_predictions -j 999 --use-conda --immediate-submit --notemp --cluster-config cluster.json --cluster 'python3 slurmSubmit.py {dependencies}'
```
#### Finally, generating figures for the vertebrate species snoRNAs is done as follows locally (of note: figures might differ slightly from those in the paper since RNAcentral data is updated frequently):
```bash
snakemake species_figures --use-conda --cores 1
```


