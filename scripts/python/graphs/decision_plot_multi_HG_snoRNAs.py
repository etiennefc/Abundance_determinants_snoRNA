#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import shap
import numpy as np

""" Create a shap decision plot containing all snoRNAs in GAS5 (for each model)"""
sno_ids_list = snakemake.params.sno_ids
expected_value_paths = snakemake.input.expected_values
shap_value_paths = snakemake.input.shap_values
model = snakemake.wildcards.models2
decision_plot_output = snakemake.output.decision_plot
colors_dict = snakemake.params.colors_dict
colors_ = [colors_dict['not_expressed'], colors_dict['not_expressed'], colors_dict['expressed']]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors_)

expected_val, shap_val = {}, {}
for i, path in enumerate(shap_value_paths):
    model_name, iteration = path.split('/')[-1].split('_shap')[0].split('_manual_')
    shap_i = pd.read_csv(path, sep='\t', index_col='gene_id_sno')
    expected_val_i = pd.read_csv(expected_value_paths[i], sep='\t')
    shap_val[iteration] = shap_i
    expected_val[iteration] = expected_val_i


def find_shap_values(shap_val_dict, expected_val_dict, sno_ids):
    """ Merge SHAP values of snoRNAs of interest (ex: all snoRNAs in same HG)
        into one list (per model) and do the same for expected values. The shap
        values are returned as a list of arrays, whereas the expected_vals are
        returned as a simple list (each element is the base value for a snoRNA
        in a given test set)."""
    shap_values, expected_vals = [None] * len(sno_ids), [None] * len(sno_ids)
    for i, sno_id in enumerate(sno_ids):
        for iteration_nb, shap_df in shap_val_dict.items():
            if sno_id in shap_df.index:  # if we find the snoRNA in that given test set
                shap_df_i = shap_df[shap_df.index == sno_id]
                shap_array_i = shap_df_i.to_numpy()
                shap_values[i] = shap_array_i
                expected_value_df = expected_val_dict[iteration_nb]
                expected_vals[i] = [expected_value_df.values[0][0]] * len(shap_df_i)
    expected_vals = [item for sublist in expected_vals for item in sublist]  # convert list of list into a simple list
    shap_values = np.concatenate(shap_values, axis=0)
    shap_values = list(shap_values[:, np.newaxis, :])

    return shap_values, expected_vals


shap_final, expected_val_final = find_shap_values(shap_val, expected_val, sno_ids_list)


# For log_reg, the model output is in log odds, so we need to convert it to probability using the logit function
col_names = shap_val['first'].columns.to_list()
col_names = [col.split('_norm')[0] for col in col_names]
if model in ['log_reg']:
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    plt.rcParams['svg.fonttype'] = 'none'
    shap.multioutput_decision_plot(expected_val_final, shap_final, 0, show=False,
                                feature_display_range=slice(-1, -50, -1), link='logit',
                                feature_names=col_names, plot_color=cmap)
    plt.savefig(decision_plot_output, bbox_inches='tight', dpi=600)
else:  # For RF and SVC, the output is already in probability, so no need to convert log odds to probability
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    plt.rcParams['svg.fonttype'] = 'none'
    shap.multioutput_decision_plot(expected_val_final, shap_final, 0, show=False,
                                feature_display_range=slice(-1, -50, -1),
                                feature_names=col_names, plot_color=cmap)
    plt.savefig(decision_plot_output, bbox_inches='tight', dpi=600)
