#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft
import numpy as np

log_reg_output = [path for path in snakemake.output.heatmap if 'log_reg' in path][0]
svc_output = [path for path in snakemake.output.heatmap if 'svc' in path][0]
rf_output = [path for path in snakemake.output.heatmap if 'rf' in path][0]

shap_values_paths = snakemake.input.shap_values
log_reg_dfs, svc_dfs, rf_dfs = [], [], []
for i, path in enumerate(shap_values_paths):
    model, iteration = path.split('/')[-1].split('_shap_')[0].rsplit('_', 1)
    df = pd.read_csv(path, sep='\t')
    df['id_model_iteration'] = df['gene_id_sno'] + f'_{model}_{iteration}'
    df = df.drop(columns='gene_id_sno')
    df = df.set_index('id_model_iteration')
    if model == 'log_reg':
        log_reg_dfs.append(df)
    elif model == 'svc':
        svc_dfs.append(df)
    elif model == 'rf':
        rf_dfs.append(df)

log_reg_concat_df = pd.concat(log_reg_dfs)
svc_concat_df = pd.concat(svc_dfs)
rf_concat_df = pd.concat(rf_dfs)


ft.heatmap_simple(log_reg_concat_df, 'plasma', 'SHAP values', log_reg_output)
ft.heatmap_simple(svc_concat_df, 'plasma', 'SHAP values', svc_output)
ft.heatmap_simple(rf_concat_df, 'plasma', 'SHAP values', rf_output)
