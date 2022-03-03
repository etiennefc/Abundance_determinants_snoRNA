#!/usr/bin/python3
import pandas as pd

all_shap_path = snakemake.input.all_shap
output_path = snakemake.output.concat_shap
sno_per_confusion_value = snakemake.input.sno_per_confusion_value
conf_val = snakemake.wildcards.confusion_value
conf_val_pair = {'FN': 'TP', 'TP': 'FN', 'FP': 'TN', 'TN': 'FP'}  # to help select only real confusion value
                                                                # (i.e. those always predicted as such across iterations and models)
conf_val_df = pd.read_csv([path for path in sno_per_confusion_value if conf_val in path][0], sep='\t')
conf_val_pair_df = pd.read_csv([path for path in sno_per_confusion_value if conf_val_pair[conf_val] in path][0], sep='\t')

# Select only real confusion_value (ex: FN) (those always predicted as such across models and iterations)
real_conf_val = list(set(conf_val_df.gene_id_sno.to_list()) - set(conf_val_pair_df.gene_id_sno.to_list()))

dfs = []
for path in all_shap_path:
    iteration = path.split('/')[-1].split('_shap')[0].split('_')[-1]
    df = pd.read_csv(path, sep='\t')
    df['iteration'] = iteration
    df = df[df['gene_id_sno'].isin(real_conf_val)]
    dfs.append(df)

concat_df = pd.concat(dfs)
concat_df.to_csv(output_path, index=False, sep='\t')
