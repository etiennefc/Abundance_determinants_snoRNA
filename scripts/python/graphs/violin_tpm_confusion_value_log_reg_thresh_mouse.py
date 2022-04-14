#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft
import numpy as np

""" Violin plot of TPM values for all confusion_value"""

confusion_value_df = pd.read_csv(snakemake.input.sno_per_confusion_value[0], sep='\t')


fn = confusion_value_df[confusion_value_df['confusion_matrix_val_log_reg_thresh'] == 'FN']
tp = confusion_value_df[confusion_value_df['confusion_matrix_val_log_reg_thresh'] == 'TP']
tn = confusion_value_df[confusion_value_df['confusion_matrix_val_log_reg_thresh'] == 'TN']
fp = confusion_value_df[confusion_value_df['confusion_matrix_val_log_reg_thresh'] == 'FP']



tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
color_dict = snakemake.params.color_dict

# Create avg_tpm and log2_avg_tpm columns in tpm_df
tpm_df['avg_tpm'] = tpm_df.filter(regex='^[A-Za-z].*_[123]$').mean(axis=1)
tpm_df['avg_tpm'] = tpm_df['avg_tpm'] + 0.0001  # add pseudocount to be able to compute log afterwards
tpm_df['gene_id_sno'] = tpm_df['gene_id']
tpm_df['log2_avg_tpm'] = np.log2(tpm_df['avg_tpm'])
tpm_df = tpm_df[['gene_id_sno', 'log2_avg_tpm', 'avg_tpm']]

# Merge each confusion value df to tpm_df
fn = fn.merge(tpm_df, how='left', on='gene_id_sno')
tp = tp.merge(tpm_df, how='left', on='gene_id_sno')
tn = tn.merge(tpm_df, how='left', on='gene_id_sno')
fp = fp.merge(tpm_df, how='left', on='gene_id_sno')

# Create the violin plot
concat_df = pd.concat([tn, tp, fn, fp])

ft.violin_wo_swarm(concat_df, 'confusion_matrix_val_log_reg_thresh', 'log2_avg_tpm', None, 'Confusion value', 'log2(average TPM)', '',
            color_dict, snakemake.output.violin)
