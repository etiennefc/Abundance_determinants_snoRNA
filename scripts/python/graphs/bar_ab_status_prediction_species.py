#!/usr/bin/python3
import pandas as pd
import functions as ft

""" Stacked bar chart of all predicted expressed/not_expressed snoRNAs per species"""
human_df = pd.read_csv(snakemake.input.human_labels, sep='\t')
mouse_df = pd.read_csv(snakemake.input.mouse_labels, sep='\t')
paths = snakemake.input.dfs
dfs = []
for i, path in enumerate(paths):
    species_name = path.split('/')[-1].split('_predicted_label')[0]
    df = pd.read_csv(path, sep='\t')
    df['species_name'] = species_name
    df = df[['predicted_label', 'species_name']]
    dfs.append(df)

# Concat dfs into 1 df
concat_df = pd.concat(dfs)


# Generate a bar chart of categorical features with a hue of gene_biotype
counts_per_feature = ft.count_list_x(concat_df, 'species_name',
                    list(snakemake.params.hue_color.keys()),
                    'predicted_label')
# Add human and mouse actual labels for comparison
human_expressed = len(human_df[human_df['abundance_cutoff_2'] == 'expressed'])
human_not_expressed = len(human_df[human_df['abundance_cutoff_2'] == 'not_expressed'])
mouse_expressed = len(mouse_df[mouse_df['abundance_cutoff'] == 'expressed'])
mouse_not_expressed = len(mouse_df[mouse_df['abundance_cutoff'] == 'not_expressed'])
counts_per_feature = [[human_expressed, human_not_expressed]] + [[mouse_expressed, mouse_not_expressed]] + counts_per_feature

# Convert to percent
percent = ft.percent_count(counts_per_feature)


# Get the total number of snoRNAs (for which we found snoRNA type) per species
total_nb_sno = str([sum(l) for l in counts_per_feature])

ft.stacked_bar2(percent, ['homo_sapiens', 'mus_musculus']+sorted(list(concat_df['species_name'].unique())),
                list(snakemake.params.hue_color.keys()), '', 'Species',
                'Proportion of snoRNAs (%)', snakemake.params.hue_color, total_nb_sno,
                snakemake.output.bar)
