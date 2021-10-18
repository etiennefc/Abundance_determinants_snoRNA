#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
import matplotlib.pyplot as plt
import shap
import numpy as np
import functions as ft

X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id_sno')
X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')
y_test = pd.read_csv(snakemake.input.y_test, sep='\t')
y_test.index = X_test.index
col_color_dict = snakemake.params.labels_dict
col_color_dict[0] = col_color_dict.pop('not_expressed')
col_color_dict[1] = col_color_dict.pop('expressed')

# Unpickle and thus instantiate the model represented by the 'models2' wildcard
# Instantiate the explainer using the X_train as backgorund data and X_test to generate shap global values
if snakemake.wildcards.models2 == "log_reg":
    model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    explainer = shap.LinearExplainer(model, shap.sample(X_train, 100, random_state=42))  # reduce number of background sample to 100
    shap_values = explainer.shap_values(X_test)

    df = pd.DataFrame(shap_values.T, index=X_test.columns, columns=X_test.index)

    # Data for heatmap column colorbar (with y_test)
    predicted_label = pd.DataFrame(model.predict(X_test))
    predicted_label.index = X_test.index
    predicted_label.columns = ['predicted_label']

    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    ft.heatmap(df, col_color_dict, y_test['label'], predicted_label['predicted_label'],
                'plasma', 'SHAP value', snakemake.output.heatmap)

else:
    model2 = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    explainer2 = shap.KernelExplainer(model2.predict, shap.sample(X_train, 100, random_state=42)) # reduce number of background sample to 100
    shap_values2 = explainer2.shap_values(X_test)

    df = pd.DataFrame(shap_values2.T, index=X_test.columns, columns=X_test.index)
    predicted_label = pd.DataFrame(model2.predict(X_test))
    predicted_label.index = X_test.index
    predicted_label.columns = ['predicted_label']

    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    ft.heatmap(df, col_color_dict, y_test['label'], predicted_label['predicted_label'],
                'plasma', 'SHAP value', snakemake.output.heatmap)
