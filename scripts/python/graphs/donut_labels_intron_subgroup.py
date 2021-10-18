#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv(snakemake.input.df, sep='\t')

# Select intronic snoRNAs and create intron subgroup (small vs long intron)
df = df[df['host_biotype2'] != 'intergenic']
df.loc[df['intron_length'] < 5000, 'intron_subgroup'] = 'small_intron'
df.loc[df['intron_length'] >= 5000, 'intron_subgroup'] = 'long_intron'

# Create a donut chart of the abundance status of snoRNAs (inner donut)
# The outer donut shows the intron subgroup of snoRNAs
count_datasets = []
count_attributes = []
for sub in snakemake.params.intron_subgroup_colors.keys():  # Iterate through abundance statuses
    temp_df = df[df['intron_subgroup'] == sub]
    count_datasets.append(len(temp_df))
    attributes_dict = {}
    for status in snakemake.params.label_colors.keys():
        attributes_dict[status] = len(temp_df[temp_df['abundance_cutoff_2'] == status])
    sorted_dict = {k: v for k, v in sorted(attributes_dict.items())}  # Sort dictionary alphabetically as in the config file
    print(sub, sorted_dict)
    for val in sorted_dict.values():
        count_attributes.append(val)

counts = [count_datasets, count_attributes]

# Set inner_labels as a list of empty strings, and labels as outer and inner_labels
inner_labels = [None] * len(snakemake.params.intron_subgroup_colors.keys()) * len(snakemake.params.label_colors.keys())
labels = [list(snakemake.params.intron_subgroup_colors.keys()), inner_labels]

# Set inner colors as a repeated list of colors for each part of the inner donut (ex: same 2 inner colors repeated for each outer donut part)
inner_colors = list(snakemake.params.label_colors.values()) * len(snakemake.params.intron_subgroup_colors.keys())

# Set colors as outer and inner_colors
colors = [list(snakemake.params.intron_subgroup_colors.values()), inner_colors]

ft.donut_2(counts, labels, colors, '', list(snakemake.params.label_colors.keys()), list(snakemake.params.label_colors.values()), snakemake.output.donut)
