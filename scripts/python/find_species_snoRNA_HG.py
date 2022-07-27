#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp
import re

species = snakemake.wildcards.species
sno_bed_path = snakemake.input.sno_bed
gtf_bed_path = snakemake.input.formatted_gtf_bed
sno_bed = BedTool(sno_bed_path)

# Intersect sno_bed with bed of genes from gtf to find if any are a snoRNA host gene
# Intersect must be on the same strand (s=True) and must fully include the snoRNA (f=1)
intersection = sno_bed.intersect(gtf_bed_path, s=True, f=1, wb=True, sorted=True).saveas(f'{species}_temp_snoRNA_HG.bed')

# Load intersect_df
cols = ['chr_sno', 'start_sno', 'end_sno', 'gene_id_sno', 'dot', 'strand_sno',
        'source_sno', 'feature_sno', 'dot2', 'attributes_sno', 'chr_host',
        'start_host', 'end_host', 'host_id', 'dot3', 'strand_host', 'source_host',
        'feature_host', 'dot4', 'attributes_host']
intersect_df = pd.read_csv(f'{species}_temp_snoRNA_HG.bed', sep='\t', names=cols, header=None)


# Retrieve host name and biotype from the attribtes_host column
attribute_dict = dict(zip(intersect_df.host_id, intersect_df.attributes_host.str.split('; ')))
host_name, host_biotype = {}, {}
for k, v in attribute_dict.items():
    if all('gene_name' not in att for att in v):  # if attribute 'gene_name' is missing in all attributes
        gene_id = [attribute for attribute in v if 'gene_id' in attribute][0]
        fake_gene_name = gene_id.replace('gene_id', 'gene_name') # create a fake gene name which will be the gene_id
        v = v + [fake_gene_name]
    name = [item.replace('"', '').replace(';', '').split(' ', maxsplit=1)[1] for item in v if 'gene_name' in item][0]
    biotype = [item.replace('"', '').replace(';', '').split(' ', maxsplit=1)[1] for item in v if 'gene_biotype' in item][0]
    host_name[k] = name
    host_biotype[k] = biotype

intersect_df['host_biotype'] = intersect_df['host_id'].map(host_biotype)
intersect_df['host_name'] = intersect_df['host_id'].map(host_name)

# Keep only relevant columns
intersect_df = intersect_df.drop(columns=['dot', 'source_sno', 'feature_sno', 'dot2',
                                'attributes_sno', 'chr_host', 'dot3', 'strand_host',
                                'source_host', 'feature_host', 'dot4', 'attributes_host'])

dfs = []
for i, group in intersect_df.groupby('gene_id_sno'):
    group['host_length'] = group['end_host'] - group['start_host'] + 1
    group = group.reset_index(drop=True)
    cols_ = group.columns
    if len(group) > 1:
        if len(group[group['host_name'].str.contains('Snhg|SNHG')]) == 1: # if only one SNHG gene is present in the potential HG, define as HG
            temp_df = group[group['host_name'].str.contains('Snhg|SNHG')]
            dfs.append(temp_df)
        elif len(group[group['host_name'].str.contains('Snhg|SNHG')]) > 1: # if multiple SNHG genes are present in the potential HG, define the shortest as HG
            temp_df = group[group['host_name'].str.contains('Snhg|SNHG')]
            temp_df = temp_df.iloc[temp_df['host_length'].idxmin()].reset_index(drop=True).to_frame()
            temp_df = temp_df.T
            temp_df.columns = cols_
            dfs.append(temp_df)
        else:
            temp_df = group.iloc[group['host_length'].idxmin()].reset_index(drop=True).to_frame()  # select the shortest potential HG
            temp_df = temp_df.T
            temp_df.columns = cols_
            dfs.append(temp_df)
    else:
        dfs.append(group)

# Concat dfs together
final_host_df = pd.concat(dfs)
final_host_df.to_csv(snakemake.output.species_snoRNA_HG, sep='\t', index=False)

sp.call(f'rm {species}_temp_snoRNA_HG.bed', shell=True)
