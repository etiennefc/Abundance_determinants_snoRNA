#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
import shap
import numpy as np

""" Compute the SHAP value of all features for all snoRNAs in each test set and each model."""

iteration_ = snakemake.wildcards.iteration
model_name = snakemake.wildcards.models2
X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id_sno')
X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')
model_path = snakemake.input.pickled_trained_model
model = pickle.load(open(model_path, 'rb'))

if model_name == 'log_reg':
    explainer = shap.LinearExplainer(model, shap.sample(X_train, 100, random_state=42))
    shap_val = explainer.shap_values(X_test)
    shap_val_df = pd.DataFrame(shap_val, index=X_test.index, columns=X_test.columns)
    shap_val_df = shap_val_df.add_suffix('_SHAP')
    shap_val_df = shap_val_df.reset_index()
    base_value = explainer.expected_value  # this is the base value where the decision plot starts (average of all X_train log odds)
    base_value_df = pd.DataFrame([base_value], columns=[f'expected_value_{model_name}_{iteration_}'])
    shap_val_df.to_csv(snakemake.output.shap, sep='\t', index=False)
    base_value_df.to_csv(snakemake.output.expected_value, sep='\t', index=False)
else:
    explainer = shap.KernelExplainer(model.predict, shap.sample(X_train, 100, random_state=42))
    shap_val = explainer.shap_values(X_test)
    shap_val_df = pd.DataFrame(shap_val, index=X_test.index, columns=X_test.columns)
    shap_val_df = shap_val_df.add_suffix('_SHAP')
    shap_val_df = shap_val_df.reset_index()
    base_value = explainer.expected_value  # this is the base value where the decision plot starts (average of all X_train log odds)
    base_value_df = pd.DataFrame([base_value], columns=[f'expected_value_{model_name}_{iteration_}'])
    shap_val_df.to_csv(snakemake.output.shap, sep='\t', index=False)
    base_value_df.to_csv(snakemake.output.expected_value, sep='\t', index=False)
