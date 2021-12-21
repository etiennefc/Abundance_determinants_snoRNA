#!/usr/bin/python3
import pandas as pd
import shap
import numpy as np
import re
""" Create a dataframe containing the rank of importance for each feature
    and per model per iteration."""

X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id_sno')
X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')
iteration = snakemake.wildcards.iteration

log_reg_model = snakemake.input.log_reg
svc_model, rf_model = snakemake.input.svc, snakemake.input.rf

# Instantiate log_reg model and get the SHAP ranking of its features (highest absolute value of SHAP is considered as rank 1 most predictive feature)
log_reg = pickle.load(open(log_reg_model[0], 'rb'))  # log_reg_model is a snakemake Namedlist of one item, so we need to select the first item through indexing
explainer_log_reg = shap.LinearExplainer(log_reg, shap.sample(X_train, 100, random_state=42))  # reduce number of background sample to 100
shap_values_log_reg = explainer_log_reg.shap_values(X_test)
vals_log_reg = np.abs(shap_values_log_reg).mean(0)  # mean SHAP value across all examples in X_test for each feature
feature_importance_log_reg = pd.DataFrame(list(zip(X_train.columns, vals_log_reg)), columns=['feature', 'feature_importance'])
feature_importance_log_reg.sort_values(by=['feature_importance'], ascending=False , inplace=True)
feature_importance_log_reg['feature_rank'] = feature_importance_log_reg.reset_index().index + 1  # Create a rank column for feature importance rank
feature_importance_log_reg['model'] = f'log_reg_{iteration}'
feature_importance_log_reg = feature_importance_log_reg.drop('feature_importance', axis=1)

# Instantiate all other 2 models (svc and rf) and get their SHAP ranking of all features
dfs = [feature_importance_log_reg]
for i, mod in enumerate([svc_model, rf_model]):
    model = pickle.load(open(mod[0], 'rb'))  # mod is a snakemake Namedlist of one item, so we need to select the first item through indexing
    model_substring = re.search("results/trained_models/(.*)_trained_scale_.*sav", mod[0]).group(1)  # find the model name within the pickled model name
    explainer = shap.KernelExplainer(model.predict, shap.sample(X_train, 100, random_state=42))  # reduce number of background sample to 100
    shap_values = explainer.shap_values(X_test)
    vals = np.abs(shap_values).mean(0)  # mean SHAP value across all examples in X_test for each feature
    feature_importance = pd.DataFrame(list(zip(X_train.columns, vals)), columns=['feature', 'feature_importance'])
    feature_importance.sort_values(by=['feature_importance'], ascending=False , inplace=True)
    feature_importance['feature_rank'] = feature_importance.reset_index().index + 1  # Create a rank column for feature importance rank
    feature_importance['model'] = f'{model_substring}_{iteration}'
    feature_importance = feature_importance.drop('feature_importance', axis=1)
    dfs.append(feature_importance)

# Concat all dfs into one df
df_final = pd.concat(dfs)
df_final.to_csv(snakemake.output.rank_features_df, sep='\t', index=False)
