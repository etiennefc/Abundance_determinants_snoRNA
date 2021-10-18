#!/usr/bin/python3
import pandas as pd


""" Fill NaN in feature df in numerical columns with -5 instead. -5 was chosen
    arbitrarily so that this negative value should not interfere with all other
    positive values. Then, apply mean normalization (feature scaling) to
    numerical features columns."""

df = pd.read_csv(snakemake.input.feature_df, sep='\t')

# Fill NaN with -5
df = df.fillna(-5)

# Scale feature values using mean normalization for numerical value columns
# with high standard deviation
df_num = df.select_dtypes(include=['int64', 'float64'])
num_cols = list(df_num.columns)
for i, col in enumerate(num_cols):
    mean = df[col].mean()
    std = df[col].std()
    df[col+'_norm'] = (df[col] - mean) / std

df = df.drop(num_cols, axis=1)
df.to_csv(snakemake.output.scaled_feature_df, index=False, sep='\t')
