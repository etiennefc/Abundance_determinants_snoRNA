#!/usr/bin/python3
import pandas as pd
import functions as ft

""" Create a bar plot per confusion value comparison (e.g. FP vs TP) for all
    categorical features in the top 10 predictive features. Each comparison counts
    only one time a snoRNA (ex: it considers a TP snoRNA once even if it is
    predicted multiple time as a TP across iterations)."""

color_dict = snakemake.params.color_dict
output = snakemake.output.bar
categorical_feature = snakemake.wildcards.top_10_categorical_features
comparison = snakemake.wildcards.comparison_confusion_val
feature_df = pd.read_csv(snakemake.input.feature_df, sep='\t')
feature_df = feature_df[['gene_id_sno', categorical_feature]]

# Get the list of all snoRNAs of a given confusion_value inside dict
sno_per_confusion_value_paths = snakemake.input.sno_per_confusion_value
sno_per_confusion_value = {}
for path in sno_per_confusion_value_paths:
    confusion_value = path.split('/')[-1]
    confusion_value = confusion_value.split('_')[0]
    df = pd.read_csv(path, sep='\t')
    sno_list = df['gene_id_sno'].to_list()
    sno_per_confusion_value[confusion_value] = sno_list


# Get the snoRNA feature value for each confusion value in the comparison
confusion_val1, confusion_val2_3 = comparison.split('_vs_')
confusion_val2, confusion_val3 = confusion_val2_3.split('_')
df1 = feature_df[feature_df['gene_id_sno'].isin(sno_per_confusion_value[confusion_val1])]
df2 = feature_df[feature_df['gene_id_sno'].isin(sno_per_confusion_value[confusion_val2])]
df3 = feature_df[feature_df['gene_id_sno'].isin(sno_per_confusion_value[confusion_val3])]

# Count feature values to create the bar chart in the order (TN, FP/FN, TP)
count_list = []
for conf_val_df in [df2, df1, df3]:
    true_condition = conf_val_df[conf_val_df[categorical_feature] == 1.0]  # ex: "host is expressed" is true
    false_condition = conf_val_df[conf_val_df[categorical_feature] == 0.0]  # ex: "host is expressed" is false
    temp = [len(true_condition), len(false_condition)]
    count_list.append(temp)

percent_count = ft.percent_count(count_list)

# Create bar plot
colors = [color_dict['True'], color_dict['False']]
ft.stacked_bar(percent_count, [confusion_val2, confusion_val1, confusion_val3],
                ["True", "False"], f'{categorical_feature}',
                'Confusion value snoRNA group', 'Proportion of snoRNAs (%)', colors, output)
