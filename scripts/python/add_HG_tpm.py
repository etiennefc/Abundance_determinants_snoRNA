#!/usr/bin/python3
import pandas as pd

"""Add host gene (HG) abundance values in each tissue to each snoRNA in a dataframe using a reference table"""

ref_table_HG = pd.read_csv(snakemake.input.ref_HG_table)
tpm_df = pd.read_csv(snakemake.input.tpm_biotype_df)
df = tpm_df[tpm_df['gene_biotype2'] == 'snoRNA']

# Merge host id, name, biotype, start, end and also sno start and end to existing dataframe
df = df.merge(ref_table_HG[['sno_id', 'sno_start', 'sno_end', 'host_id', 'host_name', 'host_biotype', 'host_start', 'host_end']],
              how='left', left_on='gene_id', right_on='sno_id')

# Create dictionary of dictionaries from tpm_df (primary key: tissue and primary value: abundance (TPM) (secondary value) for each HG id (secondary key))
tpm_df = tpm_df[tpm_df['gene_id'].isin(ref_table_HG['host_id'])].set_index('gene_id').iloc[:, 3:25].add_suffix('_host')  # select tpm col and add '_host' to col titles
tpm_dict = tpm_df.to_dict()

# Add new columns of HG tpm based on the host_id col in df and the tpm_dict
for key, val in tpm_dict.items():
    df[key] = df['host_id'].map(val)

df.to_csv(snakemake.output.sno_HG_df, index=False)
