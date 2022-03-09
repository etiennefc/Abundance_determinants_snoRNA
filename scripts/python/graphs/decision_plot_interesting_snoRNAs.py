#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import shap

""" Create a shap decision plot for SNORA77B for all models"""
output_path = snakemake.output.decision_plot
colors_dict = snakemake.params.colors_dict
colors_ = [colors_dict['not_expressed'], colors_dict['not_expressed'], colors_dict['expressed']]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors_)
sno_id = snakemake.wildcards.interesting_sno_ids
X_test_paths = snakemake.input.X_test
expected_value_paths = snakemake.input.expected_values
shap_value_paths = snakemake.input.shap_values

# Select expected_val, shap_val and feature value for the given snoRNA
expected_val, shap_val, X_test = [], [], []
for i, path in enumerate(shap_value_paths):
    model_name, iteration = path.split('/')[-1].split('_shap')[0].split('_manual_')
    shap_i = pd.read_csv(path, sep='\t', index_col='gene_id_sno')
    expected_val_i = pd.read_csv(expected_value_paths[i], sep='\t')
    X_test_i = pd.read_csv(X_test_paths[i], sep='\t', index_col='gene_id_sno')
    if sno_id in shap_i.index:
        X_test_i = X_test_i[X_test_i.index == sno_id]
        shap_i = shap_i[shap_i.index == sno_id]
        shap_i = shap_i.to_numpy()
        shap_val.append(shap_i)
        expected_val.append(expected_val_i)
        X_test.append(X_test_i)

temp_df = pd.read_csv(shap_value_paths[0], sep='\t', index_col='gene_id_sno')
col_names = temp_df.columns.to_list()
col_names = [col.split('_norm')[0] for col in col_names]

# For log_reg, the model output is in log odds, so we need to convert it to probability using the logit function
if snakemake.wildcards.models2 == "log_reg":
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    shap.decision_plot(expected_val[0].values[0][0], shap_val[0], X_test[0],
                    feature_names=col_names, show=False, feature_display_range=slice(-1, -50, -1), link='logit', plot_color=cmap)
    plt.savefig(output_path, bbox_inches='tight', dpi=600)

else: # For RF and SVC, the output is already in probability, so no need to convert log odds to probability
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    shap.decision_plot(expected_val[0].values[0][0], shap_val[0], X_test[0],
                    feature_names=col_names, show=False, feature_display_range=slice(-1, -50, -1), plot_color=cmap)
    plt.savefig(output_path, bbox_inches='tight', dpi=600)
