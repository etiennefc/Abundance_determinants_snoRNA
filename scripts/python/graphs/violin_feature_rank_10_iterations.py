#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft

df = pd.read_csv(snakemake.input.rank_features_df, sep='\t')
model_colors_dict = snakemake.params.model_colors
# Order violin plots by increasing median value of feature_ranks
ordered_violin = df.groupby('feature')['feature_rank'].median().sort_values(ascending=True).index

# Remove the iteration value from the model names (ex: log_reg instead of log_reg_first)

df['model'] = df['model'].str.rsplit('_', 1, expand=True)  # max split of 1 from the right to retain log_reg not just log

# Create the connected scatter plot
ft.violin(df, 'feature', 'feature_rank', None, 'model', 'Features', 'Predictive rank across \nmodels and iterations', '',
            ['lightgrey'], model_colors_dict, snakemake.output.violin, order=ordered_violin)
