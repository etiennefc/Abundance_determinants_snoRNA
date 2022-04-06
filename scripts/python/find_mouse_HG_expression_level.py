#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Define an abundance cutoff for host genes that will be used as a feature
    (>1 TPM in at least one average condition). The samples used to quantify HG
    abundance are from Shen et al. Nature 2012. They are mouse RNA-Seq
    experiments on stem cells, embryos and adult tissues. We select only the
    mESC and embryonic brain to use as the samples to determine the
    abundance_cutoff_host threshold. These samples were processed using the
    Recount3 pipeline. """

host_df = pd.read_csv(snakemake.input.host_df, sep='\t')
sno_tpm_df = pd.read_csv(snakemake.input.sno_tpm_df, sep='\t')
recount_tpm_df = pd.read_csv(snakemake.input.recount_tpm_df)
srr_id_dict = snakemake.params.srr_id_conversion

# Add host id column to sno_tpm_df
host_dict = dict(zip(host_df.gene_id_sno, host_df.host_id))
sno_tpm_df['host_id'] = sno_tpm_df['gene_id'].map(host_dict)

# Recount tpm df processing
recount_tpm_df[['gene_id', 'version']] = recount_tpm_df['Unnamed: 0'].str.split('.', expand=True)
recount_tpm_df = recount_tpm_df.drop(['version', 'Unnamed: 0'], axis=1)
recount_tpm_df = recount_tpm_df.rename(columns=srr_id_dict)

# Select only host genes from recount total tpm_df
recount_tpm_df = recount_tpm_df[recount_tpm_df['gene_id'].isin(list(host_df['host_id']))].reset_index(drop=True)

# Select only mESC (embryonic stem cells) and brain of mouse embryo at age E14.5
# These samples are comparable to the TGIRT-Seq samples quantifying snoRNAs in mESC and in mESC treated with retinoic acid (differentiate into neurons)
recount_tpm_df = recount_tpm_df[['gene_id', 'mESC_1', 'mESC_2', 'E14_5_brain_1', 'E14_5_brain_2']]

# If host gene is expressed >1 TPM in at least one average condition (average of the duplicates), it is expressed, else not expressed (or no host gene at all)
hg_abundance = recount_tpm_df.filter(regex='_[123]$').reset_index(drop=True)  # tpm column must end with '_1, _2 or _3'
cols_hg = list(hg_abundance.columns)
duplicates_hg = [cols_hg[n:n+2] for n in range(0, len(cols_hg), 3)]

recount_tpm_df['abundance_cutoff_host'] = ''
for i in range(0, len(recount_tpm_df)):
    row = hg_abundance.iloc[i]
    for j, duplicate in enumerate(duplicates_hg):
        if (recount_tpm_df.loc[i, duplicate].mean() > 1) == True:
            recount_tpm_df.loc[i, 'abundance_cutoff_host'] = 'host_expressed'
            break
recount_tpm_df['abundance_cutoff_host'] = recount_tpm_df['abundance_cutoff_host'].replace('', 'host_not_expressed')

# if HG was not quantified in recount tpm_df, consider it as not expressed
no_tpm_HG = [recount_tpm_df]
for host_id in host_df.host_id:
    if host_id in list(recount_tpm_df.gene_id):
        continue
    else:
        vals = [[host_id, 0, 0, 0, 0, 'host_not_expressed']]
        temp_df = pd.DataFrame(vals, columns=['gene_id', 'mESC_1', 'mESC_2', 'E14_5_brain_1', 'E14_5_brain_2', 'abundance_cutoff_host'])
        no_tpm_HG.append(temp_df)
recount_tpm_df = pd.concat(no_tpm_HG)
recount_tpm_df = recount_tpm_df.rename(columns={'gene_id': 'host_id'})

# Merge abundance_cutoff_host information to sno_tpm df and fill NA for intergenic snoRNAs
final_df = sno_tpm_df.merge(recount_tpm_df[['host_id', 'abundance_cutoff_host']], how='left', on='host_id')
final_df['abundance_cutoff_host'] = final_df['abundance_cutoff_host'].fillna('intergenic')

final_df.to_csv(snakemake.output.HG_abundance_df, index=False, sep='\t')
