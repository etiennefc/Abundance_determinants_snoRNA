#!/usr/bin/python3
import pandas as pd
import functions as ft
import numpy as np

confusion_val = snakemake.wildcards.confusion_value
confusion_value_sno_df = pd.read_csv(snakemake.input.confusion_value_df, sep='\t')
confusion_value_sno = list(confusion_value_sno_df.gene_id_sno)
shap_value_paths = snakemake.input.shap_values
shap_dfs = []

for i, path in enumerate(shap_value_paths):
    model, iteration = path.split('/')[-1].split('_shap_values')[0].rsplit('_', 1)
    df = pd.read_csv(path, sep='\t')
    df['model'] = model
    df['iteration'] = iteration
    shap_dfs.append(df)

concat_shap = pd.concat(shap_dfs)

# Keep only snoRNAs that are part of the given confusion value
concat_shap = concat_shap[concat_shap['gene_id_sno'].isin(confusion_value_sno)]

# Convert all SHAP values to absolute values of these SHAP values
concat_shap.update(concat_shap.select_dtypes(include=[np.number]).abs())

# Get the top 1, 2 and 3 feature (i.e. highest abs(SHAP value) across features) for each snoRNA
concat_shap['top_1'] = concat_shap.filter(regex="_SHAP$").apply(lambda row: row[row == row.nlargest(1).values[-1]].index[0], axis=1)
concat_shap['top_2'] = concat_shap.filter(regex="_SHAP$").apply(lambda row: row[row == row.nlargest(2).values[-1]].index[0], axis=1)
concat_shap['top_3'] = concat_shap.filter(regex="_SHAP$").apply(lambda row: row[row == row.nlargest(3).values[-1]].index[0], axis=1)

# Drop snoRNAs that have the same model and top1/2/3 (so it does not bias the global portrait
# if a snoRNA is present in multiple iterations and always predicted based on the same features)
concat_shap = concat_shap.drop_duplicates(subset=['top_1', 'top_2', 'top_3', 'model'])
concat_shap.to_csv(snakemake.output.df, sep='\t', index=False)
print(concat_shap)
len_df = len(concat_shap)

# Create a list of list containing the relative number of times a feature is classified as a top 1, 2 or 3 feature
'''
relative_number_tops = []
for feature in concat_shap.filter(regex="_SHAP$").columns:
    print(feature)
    temp = []
    for top in ['top_1', 'top_2', 'top_3']:
        number_in_top_x = len(concat_shap[concat_shap[top] == feature])
        relative_number = (number_in_top_x / len_df) * 100
        temp.append(relative_number)
    relative_number_tops.append(temp)
print(relative_number_tops)
'''
relative_number_tops = []
for top in ['top_1', 'top_2', 'top_3']:
    temp = []
    for feature in concat_shap.filter(regex="_SHAP$").columns:
        number_in_top_x = len(concat_shap[concat_shap[top] == feature])
        relative_number = (number_in_top_x / len_df) * 100
        temp.append(relative_number)
    relative_number_tops.append(temp)
print(relative_number_tops)


# Create grouped bar chart
xticklabels = [feat.split('_norm_')[0] for feat in concat_shap.filter(regex="_SHAP$").columns]
ft.bar_from_lst_of_lst(relative_number_tops, [0.15, 0.3, 0.45], ['blue', 'green', 'pink'], 0.2, xticklabels, 'Features',
        f'Proportion of all predicted {confusion_val} snoRNAs',
        ['1st most predictive feature', '2nd most predictive feature', '3rd most predictive feature'], snakemake.output.bar)
