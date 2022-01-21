#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft

df = pd.read_csv(snakemake.input.rank_features_df, sep='\t')
model_colors_dict = snakemake.params.model_colors


# Find the range of each distribution (max rank - min rank per feature) and its median
feature_distribution = {}
for i, group in enumerate(df.groupby('feature')['feature_rank']):
    feature_name = group[0]
    range_ = group[1].max() - group[1].min()
    median_ = group[1].median()
    feature_distribution[feature_name] = [median_, range_]

# Order violin plots by increasing median value of feature_ranks and by range as second sort if two features have the same median
feature_distribution_df = pd.DataFrame.from_dict(feature_distribution, columns = ['median', 'range'], orient='index')
ordered_violin = feature_distribution_df.sort_values(by=['median', 'range'], ascending=[True, True]).index

# Remove the iteration value from the model names (ex: log_reg instead of log_reg_first)
df['model'] = df['model'].str.rsplit('_', 1, expand=True)  # max split of 1 from the right to retain log_reg not just log

# Create the connected scatter plot
ft.violin(df, 'feature', 'feature_rank', None, 'model', 'Features', 'Predictive rank across \nmodels and iterations', '',
            ['lightgrey'], model_colors_dict, snakemake.output.violin, order=ordered_violin)
