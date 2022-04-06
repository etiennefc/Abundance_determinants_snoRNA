#!/usr/bin/python3
import pandas as pd

""" Get the host gene biotype for all mouse snoRNA genes (protein_coding, non_coding or intergenic)."""


host_info_df = pd.read_csv(snakemake.input.host_info_df, sep='\t') 
gtf = pd.read_csv(snakemake.input.gtf, sep='\t', names=['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'feature', 'dot2', 'attributes'])
gtf = gtf[['gene_id', 'attributes']]

print(host_info_df)
print(gtf)




