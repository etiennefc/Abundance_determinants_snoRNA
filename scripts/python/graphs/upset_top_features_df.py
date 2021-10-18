#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import shap
from upsetplot import plot as upset
import numpy as np
import re
""" Generate an upset plot of the intersesction of the top 5 most predictive
    features across all four models (wo RF)"""

X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id_sno')
X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')

# Instantiate log_reg model and get its top 5 features
log_reg = pickle.load(open(snakemake.input.log_reg, 'rb'))
explainer_log_reg = shap.LinearExplainer(log_reg, shap.sample(X_train, 100, random_state=42))  # reduce number of background sample to 100
shap_values_log_reg = explainer_log_reg.shap_values(X_test)
vals_log_reg = np.abs(shap_values_log_reg).mean(0)  # mean SHAP value across all examples in X_test for each feature
feature_importance_log_reg = pd.DataFrame(list(zip(X_train.columns, vals_log_reg)), columns=['feature', 'feature_importance'])
feature_importance_log_reg.sort_values(by=['feature_importance'], ascending=False , inplace=True)
feature_importance_log_reg['feature_rank'] = feature_importance_log_reg.reset_index().index + 1  # Create a rank column for feature importance rank (1 to 5)
feature_importance_log_reg['model'] = 'log_reg'
log_reg_df = feature_importance_log_reg.head(n=5)

# Instantiate all other 3 models (knn, gbm and svc) and get their top 5 features
dfs = [log_reg_df]
for i, mod in enumerate(snakemake.input.other_model):
    model = pickle.load(open(mod, 'rb'))
    explainer = shap.KernelExplainer(model.predict, shap.sample(X_train, 100, random_state=42))  # reduce number of background sample to 100
    shap_values = explainer.shap_values(X_test)
    vals = np.abs(shap_values).mean(0)  # mean SHAP value across all examples in X_test for each feature
    feature_importance = pd.DataFrame(list(zip(X_train.columns, vals)), columns=['feature', 'feature_importance'])
    feature_importance.sort_values(by=['feature_importance'], ascending=False , inplace=True)
    feature_importance['feature_rank'] = feature_importance.reset_index().index + 1  # Create a rank column for feature importance rank (1 to 5)
    model_substring = re.search("results/trained_models/(.*)_trained_scale_after_split.sav", mod).group(1)  # find the model name within the pickled model name
    feature_importance['model'] = model_substring  # Create a model column (model name)
    df = feature_importance.head(n=5)
    dfs.append(df)

# Concat top feature dfs into one df
all_top_features = pd.concat(dfs)
all_top_features.to_csv(snakemake.output.top_features_df, sep='\t', index=False)
