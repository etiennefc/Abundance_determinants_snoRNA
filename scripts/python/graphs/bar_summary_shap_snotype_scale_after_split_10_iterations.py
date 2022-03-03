#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
import matplotlib.pyplot as plt
import shap
import numpy as np

X_train_all = snakemake.input.X_train
X_test_all = snakemake.input.X_test
model_all = snakemake.input.pickled_trained_model
model_all = [path for path in model_all if snakemake.wildcards.models3 in path]
feature_df = pd.read_csv(snakemake.input.df, sep='\t')
cd_sno = feature_df[feature_df['sno_type'] == 'C/D']['gene_id_sno'].to_list()
haca_sno = feature_df[feature_df['sno_type'] == 'H/ACA']['gene_id_sno'].to_list()


# Unpickle and thus instantiate the model represented by the 'models3' wildcard
# Instantiate the explainer using the X_train as background data and X_test_snotype to generate shap global values
if snakemake.wildcards.models3 == "log_reg":
    shap_values_all_iterations, X_test_snotype_all_iterations = [], []
    for i, df_path in enumerate(X_train_all):
        model = pickle.load(open(model_all[i], 'rb'))
        X_train = pd.read_csv(df_path, sep='\t', index_col='gene_id_sno')
        X_test = pd.read_csv(X_test_all[i], sep='\t', index_col='gene_id_sno')
        explainer = shap.LinearExplainer(model, shap.sample(X_train, 100, random_state=42))  # number of background sample = 100

        # Split between C/D and H/ACA snoRNAs
        if snakemake.wildcards.sno_type == "CD":
            X_test_snotype = X_test[X_test.index.isin(cd_sno)]
        elif snakemake.wildcards.sno_type == "HACA":
            X_test_snotype = X_test[X_test.index.isin(haca_sno)]
        X_test_snotype_all_iterations.append(X_test_snotype)

        shap_values = explainer.shap_values(X_test_snotype)
        shap_values_all_iterations.append(shap_values)
    final_shap_values = np.concatenate(shap_values_all_iterations, axis=0)
    final_X_test_snotype = pd.concat(X_test_snotype_all_iterations)  # Concat vertically all X_test_snotype dfs to infer feature value in the summary plot
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    shap.summary_plot(final_shap_values, final_X_test_snotype, plot_type="bar", show=False, max_display=50)
    plt.savefig(snakemake.output.summary_plot, bbox_inches='tight', dpi=600)

else:
    shap_values_all_iterations, X_test_snotype_all_iterations = [], []
    for i, df_path in enumerate(X_train_all):
        model = pickle.load(open(model_all[i], 'rb'))
        X_train = pd.read_csv(df_path, sep='\t', index_col='gene_id_sno')
        X_test = pd.read_csv(X_test_all[i], sep='\t', index_col='gene_id_sno')
        explainer = shap.KernelExplainer(model.predict, shap.sample(X_train, 100, random_state=42)) # reduce number of background sample to 100

        # Split between C/D and H/ACA snoRNAs
        if snakemake.wildcards.sno_type == "CD":
            X_test_snotype = X_test[X_test.index.isin(cd_sno)]
        elif snakemake.wildcards.sno_type == "HACA":
            X_test_snotype = X_test[X_test.index.isin(haca_sno)]
        X_test_snotype_all_iterations.append(X_test_snotype)

        shap_values = explainer.shap_values(X_test_snotype)
        shap_values_all_iterations.append(shap_values)
    final_shap_values = np.concatenate(shap_values_all_iterations, axis=0)
    final_X_test_snotype = pd.concat(X_test_snotype_all_iterations)  # Concat vertically all X_test_snotype dfs to infer feature value in the summary plot
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    shap.summary_plot(final_shap_values, final_X_test_snotype, plot_type="bar", show=False, max_display=50)
    plt.savefig(snakemake.output.summary_plot, bbox_inches='tight', dpi=600)
