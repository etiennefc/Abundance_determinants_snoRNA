#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Extract snoRNAs that are encoded within the same host gene but that have
    different labels (at least one snoRNA expressed (1) and one not
    expressed(0))."""

df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
# Groupby host gene id and keep only host genes that have variable snoRNA abundance status (label)
groups = df.groupby('host_id')
dfs = []
for i, group in enumerate(groups):
    grouped_df = group[1]
    ab_status_dict = coll.Counter(grouped_df['abundance_cutoff_2'])

    if len(ab_status_dict.keys()) == 2:  # if there is at least 1 expressed and 1 not expressed snoRNA
        dfs.append(grouped_df)

# Concat all dfs together
final_df = pd.concat(dfs)
final_df = final_df[['gene_id', 'gene_name', 'host_id', 'host_name', 'abundance_cutoff_2', 'abundance_cutoff_host']]
final_df.to_csv(snakemake.output.variable_sno_labels_df, index=False, sep='\t')
