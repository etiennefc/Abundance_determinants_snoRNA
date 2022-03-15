#!/usr/bin/python3
import pandas as pd
import functions as ft

""" Create a density plot to compare all confusion values for all
    numerical features in the top 10 predictive features."""

color_dict = snakemake.params.color_dict
output = snakemake.output.density
numerical_feature = snakemake.wildcards.top_10_numerical_features
feature_df = pd.read_csv(snakemake.input.feature_df, sep='\t')
feature_df = feature_df[['gene_id_sno', numerical_feature]]
confusion_value_df = pd.read_csv(snakemake.input.sno_per_confusion_value, sep='\t')

# Get df of all snoRNAs of a given confusion_value inside dict
sno_per_confusion_value = {}
for conf_val in ['TN', 'TP', 'FN', 'FP']:
    df_temp = confusion_value_df[confusion_value_df['confusion_matrix'] == conf_val]
    sno_list = df_temp['gene_id_sno'].to_list()
    df = feature_df[feature_df['gene_id_sno'].isin(sno_list)]
    sno_per_confusion_value[conf_val] = df


# Create density plot
colors = [color_dict['TN'], color_dict['TP'], color_dict['FN'], color_dict['FP']]
dfs = [sno_per_confusion_value['TN'][numerical_feature], sno_per_confusion_value['TP'][numerical_feature],
        sno_per_confusion_value['FN'][numerical_feature], sno_per_confusion_value['FP'][numerical_feature]]
ft.density_x(dfs, numerical_feature, 'Density', 'linear', '',
            colors, ['TN', 'TP', 'FN', 'FP'], output)
