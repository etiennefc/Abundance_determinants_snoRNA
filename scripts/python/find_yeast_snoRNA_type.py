#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Find snoRNA type (C/D vs H/ACA) of yeast snoRNAs """

# Loading file and dfs
tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
snoRNA_type_df = pd.read_csv(snakemake.input.yeast_mine_df, sep='\t',
                        names=['id', 'other_id', 'gene_id', 'description'])

# Convert to lowercase the first 2 letters of gene_id (snR instead of SNR) in snoRNA_type_df
snoRNA_type_df['gene_id'] = snoRNA_type_df['gene_id'].str.replace(r'SNR([a-zA-Z0-9]*)', r'snR\1', regex=True)
snoRNA_type_df['gene_id'] = snoRNA_type_df['gene_id'].str.replace('snR17A', 'snR17a')
snoRNA_type_df['gene_id'] = snoRNA_type_df['gene_id'].str.replace('snR17B', 'snR17b')

print(snoRNA_type_df)
print(tpm_df)

# Create snoRNA_type column
snoRNA_type_df.loc[snoRNA_type_df['description'].str.contains('C/D|U3'), 'snoRNA_type'] = 'C/D'
snoRNA_type_df.loc[snoRNA_type_df['description'].str.contains('H/ACA'), 'snoRNA_type'] = 'H/ACA'
snoRNA_type_df = snoRNA_type_df.dropna(subset=['snoRNA_type'])
snoRNA_type_df = snoRNA_type_df[['gene_id', 'snoRNA_type']]

# Merge with tpm_df
sno_tpm_df = snoRNA_type_df.merge(tpm_df, how='left', on='gene_id')

sno_tpm_df.to_csv(snakemake.output.snoRNA_type_df, sep='\t', index=False)
