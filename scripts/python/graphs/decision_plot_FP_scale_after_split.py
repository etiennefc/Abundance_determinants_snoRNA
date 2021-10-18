#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import shap
import numpy as np
import subprocess as sp

""" Create a shap decision plot per specific snoRNAs that are false positives in
    all 4 models."""
output_path = snakemake.params.decision_plot_FP
log = snakemake.output.shap_local_FP_log
false_positives = snakemake.params.false_positives
sp.call("mkdir -p "+output_path+" &> "+log, shell=True)

X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id_sno')
X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')

# Unpickle and thus instantiate the model represented by the 'models' wildcard
# Instantiate the explainer using the X_train as background data and X_test to generate shap local values for one snoRNA
if snakemake.wildcards.models2 == "log_reg":
    model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    explainer = shap.LinearExplainer(model, shap.sample(X_train, 100, random_state=42))  # reduce number of background sample to 100
    for sno_id in false_positives:
        shap_values = explainer.shap_values(X_test.loc[sno_id, :])  # Select one snoRNA
        plt.rcParams['svg.fonttype'] = 'none'
        fig, ax = plt.subplots(1, 1, figsize=(15, 15))
        shap.decision_plot(explainer.expected_value, shap_values,
                        X_test.loc[sno_id, :], show=False, feature_display_range=slice(-1, -50, -1), link='logit')
        plt.savefig(output_path+sno_id+"_"+snakemake.wildcards.models2+"_all_features_test_set_100_background.svg", bbox_inches='tight', dpi=600)

else:
    model2 = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    explainer2 = shap.KernelExplainer(model2.predict, shap.sample(X_train, 100, random_state=42)) # reduce number of background sample to 100
    for sno_id in false_positives:
        shap_values2 = explainer2.shap_values(X_test.loc[sno_id, :])  # Select one snoRNA
        plt.rcParams['svg.fonttype'] = 'none'
        fig, ax = plt.subplots(1, 1, figsize=(15, 15))
        shap.decision_plot(explainer2.expected_value, shap_values2, X_test.loc[sno_id, :], show=False, feature_display_range=slice(-1, -50, -1))
        plt.savefig(output_path+sno_id+"_"+snakemake.wildcards.models2+"_all_features_test_set_100_background.svg", bbox_inches='tight', dpi=600)
