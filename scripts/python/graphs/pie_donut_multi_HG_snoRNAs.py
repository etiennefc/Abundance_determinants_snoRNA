#!/usr/bin/python3
import pandas as pd
import functions as ft

feature_df = pd.read_csv(snakemake.input.df, sep='\t')
host_df = pd.read_csv(snakemake.input.host_df)
multi_HG_different_label_snoRNAs_df = pd.read_csv(snakemake.input.multi_HG_different_label_snoRNAs, sep='\t')
merged_conf_values_paths = snakemake.params.merged_confusion_values
mono_vs_multi_HG_colors = snakemake.params.mono_vs_multi_HG_colors
multi_HG_labels_colors = snakemake.params.multi_HG_labels_colors
multi_HG_same_labels_proportion_colors = snakemake.params.multi_HG_same_labels_proportion_colors
multi_HG_diff_labels_proportion_colors = snakemake.params.multi_HG_diff_labels_proportion_colors

# Drop intergenic snoRNAs and merge host_df to all_features_df
feature_df = feature_df[feature_df['abundance_cutoff_host'] != 'intergenic']
feature_df = feature_df.merge(host_df, how='left', left_on='gene_id_sno', right_on='sno_id')

# Concat all merged confusion value dfs (one per manual iteration) into 1 df
confusion_value_dfs = []
for path in merged_conf_values_paths:
    df = pd.read_csv(path, sep='\t')
    confusion_value_dfs.append(df)
confusion_value_concat_df = pd.concat(confusion_value_dfs)

# Find the number of mono vs multi HG (respectively containing 1 vs >1 snoRNAs in the same HG)
mono_HG, multi_HG = [], []
for i, group in enumerate(feature_df.groupby('host_id')):
    grouped_df = group[1]
    host_id = group[0]
    if len(grouped_df) == 1:  # mono-HG
        mono_HG.append(host_id)
    elif len(grouped_df) > 1:  # multi-HG
        multi_HG.append(host_id)
mono_HG_nb, multi_HG_nb = len(mono_HG), len(multi_HG)
pie_counts = [mono_HG_nb, multi_HG_nb]

# Find in the multi_HG those that have snoRNAs with all the same label vs those with different labels
sno_ids_multi_HG_diff_labels = list(multi_HG_different_label_snoRNAs_df.gene_id_sno)
multi_HG_same_label = feature_df[feature_df['host_id'].isin(multi_HG)]
multi_HG_same_label = multi_HG_same_label[~multi_HG_same_label['gene_id_sno'].isin(sno_ids_multi_HG_diff_labels)]
multi_HG_same_label_nb = len(pd.unique(multi_HG_same_label['host_id']))
multi_HG_diff_label_nb = len(pd.unique(multi_HG_different_label_snoRNAs_df['host_id']))
outer_donut_counts = [multi_HG_same_label_nb, multi_HG_diff_label_nb]

# Find for multi_HG with same snoRNA labels the % of all_expressed or all_not_expressed snoRNAs
all_expressed, all_not_expressed = [], []
for i, group in enumerate(multi_HG_same_label.groupby('host_id')):
    grouped_df = group[1]
    host_id = group[0]
    if 'expressed' in list(grouped_df.abundance_cutoff_2):
        all_expressed.append(host_id)
    elif 'not_expressed' in list(grouped_df.abundance_cutoff_2):
        all_not_expressed.append(host_id)
multi_HG_all_expressed_nb = len(all_expressed)
multi_HG_all_not_expressed_nb = len(all_not_expressed)
inner_donut_same_labels_counts = [multi_HG_all_expressed_nb, multi_HG_all_not_expressed_nb]

# Find for multi_HG with different snoRNA labels the % of snoRNAs that are 50-50 expressed-not_expressed, more expressed, or more not_expressed
half_expressed_not_expressed, more_expressed, more_not_expressed = [], [], []
for i, group in enumerate(multi_HG_different_label_snoRNAs_df.groupby('host_id')):
    grouped_df = group[1]
    host_id = group[0]
    if len(grouped_df[grouped_df['abundance_cutoff_2'] == 'expressed']) == len(grouped_df[grouped_df['abundance_cutoff_2'] == 'not_expressed']):
        half_expressed_not_expressed.append(host_id)
    elif len(grouped_df[grouped_df['abundance_cutoff_2'] == 'expressed']) > len(grouped_df[grouped_df['abundance_cutoff_2'] == 'not_expressed']):
        more_expressed.append(host_id)
    elif len(grouped_df[grouped_df['abundance_cutoff_2'] == 'expressed']) < len(grouped_df[grouped_df['abundance_cutoff_2'] == 'not_expressed']):
        more_not_expressed.append(host_id)
half_expressed_not_expressed_nb = len(half_expressed_not_expressed)
more_expressed_nb = len(more_expressed)
more_not_expressed_nb = len(more_not_expressed)
inner_donut_diff_labels_counts = [half_expressed_not_expressed_nb, more_expressed_nb, more_not_expressed_nb]

# Create pie chart of mono_HG vs multi_HG
ft.pie_simple(pie_counts, mono_vs_multi_HG_colors, '', snakemake.output.pie)

## Create donut chart of multi_HG with same vs different snoRNA labels
counts_all = [outer_donut_counts, inner_donut_same_labels_counts+inner_donut_diff_labels_counts]
# Set inner_labels as a list of empty strings, and labels as outer and inner_labels
inner_labels = [None] * len(multi_HG_labels_colors.keys()) * len(multi_HG_same_labels_proportion_colors.keys()) + [None]  # + [None] is for the third inner donut label that is only present in one of the half-inner donut
labels = [list(multi_HG_labels_colors.keys()), inner_labels]
# Set inner colors as a repeated list of colors for each part of the inner donut and colors as outer and inner_colors
inner_colors = list(multi_HG_same_labels_proportion_colors.values()) + list(multi_HG_diff_labels_proportion_colors.values())
colors = [list(multi_HG_labels_colors.values()), inner_colors]
legend_labels = list(multi_HG_same_labels_proportion_colors.keys()) + list(multi_HG_diff_labels_proportion_colors.keys())
legend_colors = list(multi_HG_same_labels_proportion_colors.values()) + list(multi_HG_diff_labels_proportion_colors.values())
ft.donut_2(counts_all, labels, colors, '', legend_labels, legend_colors, snakemake.output.donut)
