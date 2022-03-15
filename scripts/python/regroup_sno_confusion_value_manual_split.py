#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Find all the snoRNAs that are either predicted as a TP, TN, FP or FN in at
    least 2 of the 3 chosen models."""

confusion_value_df_paths = snakemake.input.confusion_value_df

# Load confusion value dfs (contain confusion value for each sno in the respective test set for each model)
confusion_value_dfs = []
for path in confusion_value_df_paths:
    df = pd.read_csv(path, sep='\t')
    last_col_name = [col_name for col_name in df.columns if 'confusion_matrix_val' in col_name][0]
    df = df.rename(columns={last_col_name: last_col_name.split('_val_')[0]})
    confusion_value_dfs.append(df)


# Concat all dfs vertically
concat_df = pd.concat(confusion_value_dfs)
confusion_final_df = concat_df[['gene_id_sno', 'confusion_matrix']]

# Keep only one line per snoRNA (the chosen confusion value is the mode (i.e the most frequent value) across the 3 models prediction)
confusion_final_df = confusion_final_df.groupby('gene_id_sno')['confusion_matrix'].agg(pd.Series.mode).to_frame()
confusion_final_df = confusion_final_df.reset_index()
print(confusion_final_df)
confusion_final_df.to_csv(snakemake.output.sno_per_confusion_value, index=False, sep='\t')
