#!/usr/bin/python3
import pandas as pd
import functions as ft

""" Create a bar plot to compare confusion values for all
    categorical features in the top 10 predictive features."""

color_dict = snakemake.params.color_dict
output = snakemake.output.bar
categorical_feature = snakemake.wildcards.top_10_categorical_features
feature_df = pd.read_csv(snakemake.input.feature_df, sep='\t')
feature_df = feature_df[['gene_id_sno', categorical_feature]]
confusion_value_df = pd.read_csv(snakemake.input.sno_per_confusion_value, sep='\t')

# Get df of all snoRNAs of a given confusion_value inside dict
sno_per_confusion_value = {}
for conf_val in ['TN', 'TP', 'FN', 'FP']:
    df_temp = confusion_value_df[confusion_value_df['confusion_matrix'] == conf_val]
    sno_list = df_temp['gene_id_sno'].to_list()
    df = feature_df[feature_df['gene_id_sno'].isin(sno_list)]
    sno_per_confusion_value[conf_val] = df

dfs = [sno_per_confusion_value['TN'], sno_per_confusion_value['TP'],
        sno_per_confusion_value['FN'], sno_per_confusion_value['FP']]

# Count feature values to create the bar chart in the order (TN, FP/FN, TP)
count_list = []
for conf_val_df in dfs:
    true_condition = conf_val_df[conf_val_df[categorical_feature] == 1.0]  # ex: "host is expressed" is true
    false_condition = conf_val_df[conf_val_df[categorical_feature] == 0.0]  # ex: "host is expressed" is false
    temp = [len(true_condition), len(false_condition)]
    count_list.append(temp)

percent_count = ft.percent_count(count_list)

# Create bar plot
colors = [color_dict['True'], color_dict['False']]
ft.stacked_bar(percent_count, ['TN', 'TP', 'FN', 'FP'],
                ["True", "False"], f'{categorical_feature}',
                'Confusion value snoRNA group', 'Proportion of snoRNAs (%)', colors, output)
