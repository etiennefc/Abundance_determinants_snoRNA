#!/usr/bin/python3
import pandas as pd
import os

""" Define an abundance cutoff for host genes that will be used as a feature
    (>1 TPM in at least one average condition). The samples used to quantify HG
    abundance are from the Bgee database of normal animal samples (Bastian et
    al. 2021). """

host_df = pd.read_csv(snakemake.input.host_df, sep='\t')
sno_df = pd.read_csv(snakemake.input.sno_df, sep='\t')
tpm_df_dir = snakemake.input.tpm_df_dir

# Add host id column to sno_df
host_dict = dict(zip(host_df.gene_id_sno, host_df.host_id))
sno_df['host_id'] = sno_df['gene_id'].map(host_dict)


# Load all TPM dfs and filter to keep only the gene_id, condition and TPM value of host genes
dfs = []
host_ids = list(pd.unique(host_df.host_id))
for file in os.listdir(tpm_df_dir):
    if file.endswith('.tsv'):
        df = pd.read_csv(f'{tpm_df_dir}/{file}', sep='\t')
        df = df[['Gene ID', 'Anatomical entity name', 'TPM']]
        df.columns = ['gene_id', 'condition', 'TPM']
        df = df[df['gene_id'].isin(host_ids)]
        dfs.append(df)

# Concat all dfs and groupby gene_id and condition; return average TPM per gene_id and condition
concat_df = pd.concat(dfs)
grouped_df = concat_df.groupby(['gene_id', 'condition'])['TPM'].mean()
grouped_df = grouped_df.reset_index()

# If host gene is expressed >1 TPM in at least one average condition (average of the replicates), it is expressed,
# else not expressed
pivot_df = grouped_df.pivot_table(index='gene_id', columns='condition', values='TPM')
pivot_df.loc[(pivot_df.iloc[:, :] > 1).any(axis=1), 'abundance_cutoff_host'] = 'host_expressed'
pivot_df['abundance_cutoff_host'] = pivot_df.abundance_cutoff_host.fillna('host_not_expressed')
pivot_df = pivot_df.reset_index()

# Merge pivot_df to host_df (if host not present/quantified in Bgee datasets, consider it as host_not_expressed)
host_df = host_df[['host_id', 'host_name', 'host_biotype']]
host_merge_df = host_df.drop_duplicates(subset='host_id').merge(pivot_df, how='left', left_on='host_id', right_on='gene_id')
host_merge_df['abundance_cutoff_host'] = host_merge_df['abundance_cutoff_host'].fillna('host_not_expressed')
host_merge_df = host_merge_df.drop(columns='gene_id')

# Merge host_abundance_cutoff info to sno_df
final_df = sno_df.merge(host_merge_df, how='left', on='host_id')
final_df['abundance_cutoff_host'] = final_df['abundance_cutoff_host'].fillna('intergenic')
temp_df = final_df.pop('abundance_cutoff_host')
final_df.insert(4, temp_df.name, temp_df)  # move abundance_cutoff_host column to the fifth position in df

final_df.to_csv(snakemake.output.HG_abundance_df, index=False, sep='\t')
