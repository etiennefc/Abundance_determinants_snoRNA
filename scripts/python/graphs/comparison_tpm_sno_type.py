#!/usr/bin/python3
import pandas as pd
import functions as ft
from scipy import stats as st
import numpy as np

tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
feature_df = pd.read_csv(snakemake.input.all_features, sep='\t')

# Create average TPM (and log10) columns in tpm_df and merge that column to feature_df
tpm_df['avg_tpm'] = tpm_df.filter(regex='^[A-Z].*_[1-3]$', axis=1).mean(axis=1)
tpm_df = tpm_df[['gene_id', 'avg_tpm']]
feature_df = feature_df.merge(tpm_df, how='left', left_on='gene_id_sno', right_on='gene_id')
feature_df = feature_df.drop(['gene_id'], axis=1)
feature_df['avg_tpm_log10'] = np.log10(feature_df['avg_tpm'])

# Get expressed snoRNAs
expressed = feature_df[feature_df['abundance_cutoff_2'] == "expressed"]


# Create violin plots to compare the abundance of expressed snoRNAs per snoRNA type
ft.violin(expressed, "sno_type", "avg_tpm_log10", None, None, "Type of snoRNA",
                "Average abundance across \n tissues (log10(TPM))", "",
                snakemake.params.colors, ['black'], snakemake.output.violin)
