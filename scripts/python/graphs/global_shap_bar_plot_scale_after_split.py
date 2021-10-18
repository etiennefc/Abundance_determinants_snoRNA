#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import shap
import numpy as np

X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id_sno')
X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')

# Unpickle and thus instantiate the model represented by the 'models' wildcard
# Instantiate the explainer using the X_train as background data and X_test to generate shap local values for one snoRNA
if snakemake.wildcards.models2 == "log_reg":
    model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    explainer = shap.LinearExplainer(model, shap.sample(X_train, 100, random_state=42))  # reduce number of background sample to 100
    shap_values = explainer.shap_values(X_test)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    shap.summary_plot(shap_values, X_test, plot_type='bar', max_display=50, show=False)
    plt.savefig(snakemake.output.bar_plot, bbox_inches='tight', dpi=600)

else:
    model2 = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    explainer2 = shap.KernelExplainer(model2.predict, shap.sample(X_train, 100, random_state=42)) # reduce number of background sample to 100
    shap_values2 = explainer2.shap_values(X_test)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    shap.summary_plot(shap_values2, X_test, plot_type='bar', max_display=50, show=False)
    plt.savefig(snakemake.output.bar_plot, bbox_inches='tight', dpi=600)
