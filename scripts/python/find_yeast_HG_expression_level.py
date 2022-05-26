#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Define an abundance cutoff for host genes that will be used as a feature
    (>1 TPM in at least one average condition). The samples used to quantify HG
    abundance are 3 WT S. cerevisiae samples processed using the TGIRT-Seq
    pipeline (Fafard-Couture et al., 2021, Genome Biology)."""

host_df = pd.read_csv(snakemake.input.host_df, sep='\t')
sno_tpm_df = pd.read_csv(snakemake.input.sno_tpm_df, sep='\t')
tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')

# Add host id column to sno_tpm_df
host_dict = dict(zip(host_df.gene_id_sno, host_df.host_name))
sno_tpm_df['host_id'] = sno_tpm_df['gene_id'].map(host_dict)
# Select only host genes from total tpm_df
tpm_df = tpm_df[tpm_df['gene_id'].isin(list(host_df['host_id']))].reset_index(drop=True)

# If host gene is expressed >1 TPM in at least one average condition (average of the duplicates), it is expressed, else not expressed (or no host gene at all)
hg_abundance = tpm_df.filter(regex='_[123]$').reset_index(drop=True)  # tpm column must end with '_1, _2 or _3'
cols_hg = list(hg_abundance.columns)
triplicates_hg = [cols_hg[n:n+3] for n in range(0, len(cols_hg), 3)]

tpm_df['abundance_cutoff_host'] = ''
for i in range(0, len(tpm_df)):
    row = hg_abundance.iloc[i]
    for j, duplicate in enumerate(triplicates_hg):
        if (tpm_df.loc[i, duplicate].mean() > 1) == True:
            tpm_df.loc[i, 'abundance_cutoff_host'] = 'host_expressed'
            break
tpm_df['abundance_cutoff_host'] = tpm_df['abundance_cutoff_host'].replace('', 'host_not_expressed')

# Merge abundance_cutoff_host information to sno_tpm df and fill NA for intergenic snoRNAs
final_df = sno_tpm_df.merge(tpm_df[['gene_name', 'abundance_cutoff_host']], how='left', left_on='host_id', right_on='gene_name')
final_df['abundance_cutoff_host'] = final_df['abundance_cutoff_host'].fillna('intergenic')
final_df = final_df.drop(columns='gene_name_y')
final_df = final_df.rename(columns={'gene_name_x': 'gene_name'})
final_df.to_csv(snakemake.output.HG_abundance_df, index=False, sep='\t')
