#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
import matplotlib.pyplot as plt
import shap
import numpy as np

""" Create a clustered shap decision plot containing all snoRNAs of a given confusion value per model."""
log_reg_output, svc_output, rf_output = snakemake.output.shap_local_log_reg, snakemake.output.shap_local_svc, snakemake.output.shap_local_rf
conf_val = snakemake.wildcards.confusion_value
sno_per_confusion_value = snakemake.input.sno_per_confusion_value
conf_val_pair = {'FN': 'TP', 'TP': 'FN', 'FP': 'TN', 'TN': 'FP'}  # to help select only real confusion value
                                                                # (i.e. those always predicted as such across iterations and models)
conf_val_df = pd.read_csv([path for path in sno_per_confusion_value if conf_val in path][0], sep='\t')
conf_val_pair_df = pd.read_csv([path for path in sno_per_confusion_value if conf_val_pair[conf_val] in path][0], sep='\t')

# Load shap_val and expe dfs into dict for all models and iteration
shap_val, expected_val = {}, {}
shap_val_paths, expected_val_paths = snakemake.input.shap_values, snakemake.input.expected_value
for i, path in enumerate(shap_val_paths):
    model_name, iteration = path.split('/')[-1].split('_shap')[0].rsplit('_', maxsplit=1)
    shap_i = pd.read_csv(path, sep='\t', index_col='gene_id_sno')
    expected_val_i = pd.read_csv(expected_val_paths[i], sep='\t')
    if model_name not in shap_val.keys():
        shap_val[model_name] = {iteration: shap_i}
    else:
        shap_val[model_name][iteration] = shap_i
    if model_name not in expected_val.keys():
        expected_val[model_name] = {iteration: expected_val_i}
    else:
        expected_val[model_name][iteration] = expected_val_i



# Select only real confusion_value (ex: FN) (those always predicted as such across models and iterations)
real_conf_val = list(set(conf_val_df.gene_id_sno.to_list()) - set(conf_val_pair_df.gene_id_sno.to_list()))

# Get shap values and expected values in the desired format for the decision plot
def merge_shap_values(shap_val_dict, expected_val_dict, model_name_, real_confusion_values):
    """ Merge SHAP values of 10 iterations into one list (per model) and do the
        same for expected values. The shap values are returned as a list of arrays
        (each array corresponds to all the snoRNA of a given confusion_value in
        one test set (iteration)), whereas the expected_vals are returned as a
        simple list (each element is the base value for a snoRNA in a given test set)"""
    shap_values = [None] * len(shap_val_dict[model_name_].keys())
    expected_vals = [None] * len(shap_val_dict[model_name_].keys())
    for i, iteration_nb in enumerate(shap_val_dict[model_name_].keys()):
        df_shap_i = shap_val_dict[model_name_][iteration_nb]
        df_shap_i = df_shap_i[df_shap_i.index.isin(real_confusion_values)]
        conf_val_sno_nb = len(df_shap_i)
        shap_i_array = df_shap_i.to_numpy()
        shap_values[i] = shap_i_array
        expected_value_df = expected_val_dict[model_name_][iteration_nb]
        expected_vals[i] =  [expected_value_df.values[0][0]] * conf_val_sno_nb
    expected_vals = [item for sublist in expected_vals for item in sublist]  # convert list of list into a simple list
    shap_values = np.concatenate(shap_values, axis=0)
    shap_values = list(shap_values[:, np.newaxis, :])
    return shap_values, expected_vals

shap_log_reg, expected_val_log_reg = merge_shap_values(shap_val, expected_val, 'log_reg', real_conf_val)
shap_svc, expected_val_svc = merge_shap_values(shap_val, expected_val, 'svc', real_conf_val)
shap_rf, expected_val_rf = merge_shap_values(shap_val, expected_val, 'rf', real_conf_val)

# For log_reg, the model output is in log odds, so we need to convert it to probability using the logit function
fig, ax = plt.subplots(1, 1, figsize=(15, 15))
shap.multioutput_decision_plot(expected_val_log_reg, shap_log_reg, 0, show=False,
                            feature_display_range=slice(-1, -50, -1), link='logit', feature_order='hclust',
                            feature_names=shap_val['log_reg']['first'].columns.to_list())
plt.savefig(log_reg_output, bbox_inches='tight', dpi=600)

# For RF and SVC, the output is already in probability, so no need to convert log odds to probability
fig, ax = plt.subplots(1, 1, figsize=(15, 15))
shap.multioutput_decision_plot(expected_val_svc, shap_svc, 0, show=False,
                            feature_display_range=slice(-1, -50, -1), feature_order='hclust',
                            feature_names=shap_val['log_reg']['first'].columns.to_list())
plt.savefig(svc_output, bbox_inches='tight', dpi=600)

fig, ax = plt.subplots(1, 1, figsize=(15, 15))
shap.multioutput_decision_plot(expected_val_rf, shap_rf, 0, show=False,
                            feature_display_range=slice(-1, -50, -1), feature_order='hclust',
                            feature_names=shap_val['log_reg']['first'].columns.to_list())
plt.savefig(rf_output, bbox_inches='tight', dpi=600)
