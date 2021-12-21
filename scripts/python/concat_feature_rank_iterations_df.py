#!/usr/bin/python3
import pandas as pd

""" Concat all iterations dfs (of feature rank per model) into one df."""

dfs = []
for i, df in enumerate(snakemake.input.dfs):
    temp_df = pd.read_csv(df, sep='\t')
    dfs.append(temp_df)

# Concat all dfs into one df
df_final = pd.concat(dfs)
df_final.to_csv(snakemake.output.concat_df, sep='\t', index=False)
