#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Find all the snoRNAs that are either a TP, TN, FP or FN in all the 10
    iterations and regroup them respectively in 4 datasets."""

confusion_value = snakemake.wildcards.confusion_value
confusion_value_df_paths = snakemake.input.confusion_value_df

# Load confusion value dfs (contain confusion value for each sno in the respective test set)
confusion_value_dfs = []
for path in confusion_value_df_paths:
    df = pd.read_csv(path, sep='\t')
    last_col_name = [col_name for col_name in df.columns if 'confusion_matrix_val' in col_name][0]
    df = df.rename(columns={last_col_name: last_col_name.split('_val_')[0]})
    confusion_value_dfs.append(df)


# Concat all dfs vertically and create one df per confusion_value
concat_df = pd.concat(confusion_value_dfs)
confusion_final_df = concat_df[concat_df['confusion_matrix'] == confusion_value]
confusion_final_df = confusion_final_df[['gene_id_sno', 'confusion_matrix']]
print(confusion_final_df)
confusion_final_df.to_csv(snakemake.output.sno_per_confusion_value, index=False, sep='\t')
