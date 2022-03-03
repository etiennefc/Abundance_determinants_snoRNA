#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft
import numpy as np


shap_values_paths = snakemake.input.shap_values
dfs = []
for i, path in enumerate(shap_values_paths):
    model, iteration = path.split('/')[-1].split('_shap_')[0].rsplit('_', 1)
    df = pd.read_csv(path, sep='\t')
    df['id_model_iteration'] = df['gene_id_sno'] + f'_{model}_{iteration}'
    df = df.drop(columns='gene_id_sno')
    df = df.set_index('id_model_iteration')
    dfs.append(df)

concat_df = pd.concat(dfs)


ft.heatmap_simple(concat_df, 'plasma', 'SHAP values', snakemake.output.heatmap)
