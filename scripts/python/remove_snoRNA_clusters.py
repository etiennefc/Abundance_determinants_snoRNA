#!/usr/bin/python3
import pandas as pd

""" Remove snoRNA clusters in SNHG14 and MEG8 host genes (so mostly SNORD115,
    116 and 113, 114 snoRNAs) from the orginal dataset containing all snoRNAs. """

ref_table_HG = pd.read_csv(snakemake.input.host_gene_df)
df = pd.read_csv(snakemake.input.df, sep='\t')

# Get all snoRNAs in SNHG14 and MEG8
clusters = ref_table_HG[(ref_table_HG['host_id'] == 'ENSG00000224078') | (ref_table_HG['host_id'] == 'ENSG00000225746')]
cluster_sno = list(clusters['sno_id'])

# Remove these snoRNAs from the original dataset
df = df[~df['gene_id_sno'].isin(cluster_sno)]

df.to_csv(snakemake.output.df_wo_clusters, index=False, sep='\t')
