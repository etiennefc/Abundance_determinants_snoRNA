#!/usr/bin/python3
import pandas as pd
import numpy as np
""" Create a dataframe containing the rank of importance for each feature
    and per model per iteration."""
manual_iteration = snakemake.wildcards.manual_iteration
shap_paths = snakemake.input.shap_vals
log_reg_shap_path = [path for path in shap_paths if 'log_reg' in path][0]
svc_shap_path = [path for path in shap_paths if 'svc' in path][0]
rf_shap_path = [path for path in shap_paths if 'rf' in path][0]

shap_values_log_reg = pd.read_csv(log_reg_shap_path, sep='\t')
shap_values_log_reg = shap_values_log_reg.set_index('gene_id_sno')
log_reg_cols = list(shap_values_log_reg.columns)
log_reg_cols = [col.split('_norm_SHAP')[0] for col in log_reg_cols]
shap_values_log_reg.columns = log_reg_cols # remove the _norm_SHAP suffix
vals_log_reg = np.abs(shap_values_log_reg).mean(0)  # mean SHAP value across all examples in X_test for each feature
feature_importance_log_reg = pd.DataFrame(list(zip(shap_values_log_reg.columns, vals_log_reg)), columns=['feature', 'feature_importance'])
feature_importance_log_reg.sort_values(by=['feature_importance'], ascending=False , inplace=True)
feature_importance_log_reg['feature_rank'] = feature_importance_log_reg.reset_index().index + 1  # Create a rank column for feature importance rank
feature_importance_log_reg['model'] = f'log_reg_{manual_iteration}'
feature_importance_log_reg = feature_importance_log_reg.drop('feature_importance', axis=1)

# Get the SHAP ranking of all features for the svc and rf models
shap_values_svc = pd.read_csv(svc_shap_path, sep='\t')
shap_values_svc = shap_values_svc.set_index('gene_id_sno')
svc_cols = list(shap_values_svc.columns)
svc_cols = [col.split('_norm_SHAP')[0] for col in svc_cols]
shap_values_svc.columns = svc_cols # remove the _norm_SHAP suffix
shap_values_rf = pd.read_csv(rf_shap_path, sep='\t')
shap_values_rf = shap_values_rf.set_index('gene_id_sno')
rf_cols = list(shap_values_rf.columns)
rf_cols = [col.split('_norm_SHAP')[0] for col in rf_cols]
shap_values_rf.columns = rf_cols # remove the _norm_SHAP suffix

dfs = [feature_importance_log_reg]
paths = [svc_shap_path, rf_shap_path]
for i, shap_values in enumerate([shap_values_svc, shap_values_rf]):
    specific_path = paths[i]
    model_substring = specific_path.split('/')[-1].split('_manual')[0]  # find the model name within the path
    vals = np.abs(shap_values).mean(0)  # mean SHAP value across all examples in X_test for each feature
    feature_importance = pd.DataFrame(list(zip(shap_values.columns, vals)), columns=['feature', 'feature_importance'])
    feature_importance.sort_values(by=['feature_importance'], ascending=False , inplace=True)
    feature_importance['feature_rank'] = feature_importance.reset_index().index + 1  # Create a rank column for feature importance rank
    feature_importance['model'] = f'{model_substring}_{manual_iteration}'
    feature_importance = feature_importance.drop('feature_importance', axis=1)
    dfs.append(feature_importance)

# Concat all dfs into one df
df_final = pd.concat(dfs)
df_final.to_csv(snakemake.output.rank_features_df, sep='\t', index=False)
