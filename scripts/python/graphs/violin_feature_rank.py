#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft

df = pd.read_csv(snakemake.input.rank_features_df, sep='\t')
model_colors_dict = snakemake.params.model_colors
# Order violin plots by increasing median value of feature_ranks
ordered_violin = df.groupby('feature')['feature_rank'].median().sort_values(ascending=True).index

# Create the connected scatter plot
ft.violin(df, 'feature', 'feature_rank', None, 'model', 'Feature', 'Predictive rank \n across models', '',
            ['lightgrey'], model_colors_dict, snakemake.output.violin, order=ordered_violin)
