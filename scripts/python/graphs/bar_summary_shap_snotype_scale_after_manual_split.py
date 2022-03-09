#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import shap
import numpy as np

X_test_paths = snakemake.input.X_test
shap_iterations_paths = snakemake.input.shap_values
feature_df = pd.read_csv(snakemake.input.df, sep='\t')
cd_sno = feature_df[feature_df['sno_type'] == 'C/D']['gene_id_sno'].to_list()
haca_sno = feature_df[feature_df['sno_type'] == 'H/ACA']['gene_id_sno'].to_list()

# Load all manual split iterations dfs and select only C/D or H/ACA (shap values and feature values)
shap_values_all_iterations, X_test_snotype_all_iterations = [], []
for i, df_path in enumerate(X_test_paths):
    X_test = pd.read_csv(X_test_paths[i], sep='\t', index_col='gene_id_sno')
    shap_iteration = pd.read_csv(shap_iterations_paths[i], sep='\t', index_col='gene_id_sno')

    # Split between C/D and H/ACA snoRNAs
    if snakemake.wildcards.sno_type == "CD":
        X_test_snotype = X_test[X_test.index.isin(cd_sno)]
        shap_iteration_sno_type = shap_iteration[shap_iteration.index.isin(cd_sno)]
    elif snakemake.wildcards.sno_type == "HACA":
        X_test_snotype = X_test[X_test.index.isin(haca_sno)]
        shap_iteration_sno_type = shap_iteration[shap_iteration.index.isin(haca_sno)]

    X_test_snotype_all_iterations.append(X_test_snotype)
    shap_values_all_iterations.append(shap_iteration_sno_type)

# Concat values of all 10 iterations in a df
final_shap_values = np.concatenate(shap_values_all_iterations, axis=0)
final_X_test_snotype = pd.concat(X_test_snotype_all_iterations)  # Concat vertically all X_test_snotype dfs to infer feature value in the summary plot

# Create summary bar plot
plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(15, 15))
shap.summary_plot(final_shap_values, final_X_test_snotype, plot_type="bar", show=False, max_display=50, color="#969696")
plt.savefig(snakemake.output.summary_plot, bbox_inches='tight', dpi=600)
