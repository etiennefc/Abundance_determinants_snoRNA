#!/usr/bin/python3
import pandas as pd

feature_df = pd.read_csv(snakemake.input.all_features_df, sep='\t')
output_path = snakemake.output.real_confusion_value_df
sno_per_confusion_value = snakemake.input.sno_per_confusion_value
conf_val = snakemake.wildcards.confusion_value
conf_val_pair = {'FN': 'TP', 'TP': 'FN', 'FP': 'TN', 'TN': 'FP'}  # to help select only real confusion value
                                                                # (i.e. those always predicted as such across iterations and models)
conf_val_df = pd.read_csv([path for path in sno_per_confusion_value if conf_val in path][0], sep='\t')
conf_val_pair_df = pd.read_csv([path for path in sno_per_confusion_value if conf_val_pair[conf_val] in path][0], sep='\t')

# Select only real confusion_value (ex: FN) (those always predicted as such across models and iterations)
real_conf_val = list(set(conf_val_df.gene_id_sno.to_list()) - set(conf_val_pair_df.gene_id_sno.to_list()))

df = feature_df[feature_df['gene_id_sno'].isin(real_conf_val)]
df.to_csv(output_path, index=False, sep='\t')
