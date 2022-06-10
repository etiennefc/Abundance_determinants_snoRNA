#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv(snakemake.input.df, sep='\t')
host_biotype_df = pd.read_csv(snakemake.input.host_biotype_df, sep='\t')
host_biotype_df = host_biotype_df[['gene_id_sno', 'host_biotype']]
host_biotype_df = host_biotype_df.rename(columns={'gene_id_sno': 'gene_id'})

simplified_biotype_dict = {'lncRNA': 'non_coding', 'protein_coding': 'protein_coding', 'TEC': 'non_coding', 'unitary_pseudogene': 'non_coding', 'unprocessed_pseudogene': 'non_coding'}
host_biotype_df['host_biotype2'] = host_biotype_df['host_biotype'].map(simplified_biotype_dict)

# Drop duplicates (keep only 1 occurence)
df = df.drop_duplicates(subset=['sno_mfe', 'terminal_stem_mfe',
                                'combined_box_hamming',
                                'abundance_cutoff_host'])

# Merge host biotype df to sno df
df = df.merge(host_biotype_df, how='left', left_on='gene_id_sno', right_on='gene_id')
df['host_biotype2'] = df['host_biotype2'].fillna('intergenic')

# Create a donut chart of the abundance status of snoRNAs (outer donut)
# The inner donut shows the host biotype (intergenic, protein-coding or non-coding)

count_datasets = []
count_attributes = []
for status in snakemake.params.label_colors.keys():  # Iterate through abundance statuses
    temp_df = df[df['abundance_cutoff'] == status]
    count_datasets.append(len(temp_df))
    attributes_dict = {}
    for type in snakemake.params.host_biotype_colors.keys():
        attributes_dict[type] = len(temp_df[temp_df['host_biotype2'] == type])
    sorted_dict = {k: v for k, v in sorted(attributes_dict.items())}  # Sort dictionary alphabetically as in the config file
    print(status, sorted_dict)
    for val in sorted_dict.values():
        count_attributes.append(val)

counts = [count_datasets, count_attributes]

# Set inner_labels as a list of empty strings, and labels as outer and inner_labels
inner_labels = [None] * len(snakemake.params.label_colors.keys()) * len(snakemake.params.host_biotype_colors.keys())
labels = [list(snakemake.params.label_colors.keys()), inner_labels]

# Set inner colors as a repeated list of colors for each part of the inner donut (ex: same 2 inner colors repeated for each outer donut part)
inner_colors = list(snakemake.params.host_biotype_colors.values()) * len(snakemake.params.label_colors.keys())

# Set colors as outer and inner_colors
colors = [list(snakemake.params.label_colors.values()), inner_colors]

ft.donut_2(counts, labels, colors, '', list(snakemake.params.host_biotype_colors.keys()), list(snakemake.params.host_biotype_colors.values()), snakemake.output.donut)
