#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv(snakemake.input.df, sep='\t')
snoRNA_type_df = pd.read_csv(snakemake.input.snoRNA_type_df, sep='\t')

# Merge dfs
df = df.merge(snoRNA_type_df, how='left', left_on='gene_id_sno', right_on='gene_id')

# Create a donut chart of the predicted abundance status of snoRNAs (outer donut)
# The inner donut shows the snoRNA type (C/D or H/ACA)

count_datasets = []
count_attributes = []
for status in snakemake.params.label_colors.keys():  # Iterate through abundance statuses
    temp_df = df[df['predicted_label'] == status]
    count_datasets.append(len(temp_df))
    attributes_dict = {}
    for type in snakemake.params.sno_type_colors.keys():
        attributes_dict[type] = len(temp_df[temp_df['snoRNA_type'] == type])
    sorted_dict = {k: v for k, v in sorted(attributes_dict.items())}  # Sort dictionary alphabetically as in the config file
    print(status, sorted_dict)
    for val in sorted_dict.values():
        count_attributes.append(val)

counts = [count_datasets, count_attributes]

# Set inner_labels as a list of empty strings, and labels as outer and inner_labels
inner_labels = [None] * len(snakemake.params.label_colors.keys()) * len(snakemake.params.sno_type_colors.keys())
labels = [list(snakemake.params.label_colors.keys()), inner_labels]

# Set inner colors as a repeated list of colors for each part of the inner donut (ex: same 2 inner colors repeated for each outer donut part)
inner_colors = list(snakemake.params.sno_type_colors.values()) * len(snakemake.params.label_colors.keys())

# Set colors as outer and inner_colors
colors = [list(snakemake.params.label_colors.values()), inner_colors]

ft.donut_2(counts, labels, colors, 'SnoRNA type according to the\npredicted abundance status', list(snakemake.params.sno_type_colors.keys()), list(snakemake.params.sno_type_colors.values()), snakemake.output.donut)
