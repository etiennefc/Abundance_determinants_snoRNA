#!/usr/bin/python3
import pandas as pd
import functions as ft

color_dict = snakemake.params.color_dict
dfs_outputs = snakemake.output.dfs
colors = [color_dict['TN'], color_dict['TP'], color_dict['FN'], color_dict['FP']]
confusion_value_paths = snakemake.input.confusion_value_dfs
rbp_enrichment_df = pd.read_csv(snakemake.input.combined_rbp_score_df, sep='\t')
confusion_value_dict = {}
for path in confusion_value_paths:
    confusion_value = path.split("/")[-1].split("_w_")[0]
    print(confusion_value)
    df = pd.read_csv(path, sep='\t')
    df = df.merge(rbp_enrichment_df, how='left', on='gene_id_sno')
    output_path = [output for output in dfs_outputs if confusion_value in output][0]
    df.to_csv(output_path, sep='\t', index=False)
    confusion_value_dict[confusion_value] = df

df_list = [confusion_value_dict['TN'].combined_rbp_score_log10,
            confusion_value_dict['TP'].combined_rbp_score_log10,
            confusion_value_dict['FN'].combined_rbp_score_log10,
            confusion_value_dict['FP'].combined_rbp_score_log10]
ft.density_x(df_list, 'Combined RBP enrichment score (log10)', 'Density', 'linear', '',
                colors, ['TN', 'TP', 'FN', 'FP'], snakemake.output.density)
