#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
import matplotlib.pyplot as plt
import shap
import numpy as np

X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id_sno')
X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')
feature_df = pd.read_csv(snakemake.input.df, sep='\t')
cd_sno = feature_df[feature_df['sno_type'] == 'C/D']['gene_id_sno'].to_list()
haca_sno = feature_df[feature_df['sno_type'] == 'H/ACA']['gene_id_sno'].to_list()

# Split between C/D and H/ACA snoRNAs
if snakemake.wildcards.sno_type == "CD":
    X_test_snotype = X_test[X_test.index.isin(cd_sno)]
elif snakemake.wildcards.sno_type == "HACA":
    X_test_snotype = X_test[X_test.index.isin(haca_sno)]
print(X_test_snotype)
print(len(X_test_snotype))

# Unpickle and thus instantiate the model represented by the 'models2' wildcard
# Instantiate the explainer using the X_train as backgorund data and X_test to generate shap global values
if snakemake.wildcards.models2 == "log_reg":
    model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    #explainer = shap.LinearExplainer(model, X_train)  # Use whole X_train as background (quite longer than using the line below with subsampled background)
    explainer = shap.LinearExplainer(model, shap.sample(X_train, 100, random_state=42))  # reduce number of background sample to 100
    shap_values = explainer.shap_values(X_test_snotype)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    shap.summary_plot(shap_values, X_test_snotype, show=False, max_display=50)
    plt.savefig(snakemake.output.summary_plot, bbox_inches='tight', dpi=600)

else:
    model2 = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    #explainer2 = shap.KernelExplainer(model2.predict, X_train)  # Use whole X_train as background (quite longer than using the line below with subsampled background)
    explainer2 = shap.KernelExplainer(model2.predict, shap.sample(X_train, 100, random_state=42)) # reduce number of background sample to 100
    shap_values2 = explainer2.shap_values(X_test_snotype)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    shap.summary_plot(shap_values2, X_test_snotype, show=False, max_display=50)
    plt.savefig(snakemake.output.summary_plot, bbox_inches='tight', dpi=600)
