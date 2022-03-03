#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft
import numpy as np

""" Violin plot of TPM values for all real False negatives and real True Positives"""

confusion_value_df_paths = snakemake.input.sno_per_confusion_value
fn_path = [path for path in confusion_value_df_paths if 'FN' in path][0]
tp_path = [path for path in confusion_value_df_paths if 'TP' in path][0]
tn_path = [path for path in confusion_value_df_paths if 'TN' in path][0]
fp_path = [path for path in confusion_value_df_paths if 'FP' in path][0]

fn = pd.read_csv(fn_path, sep='\t')
tp = pd.read_csv(tp_path, sep='\t')
tn = pd.read_csv(tn_path, sep='\t')
fp = pd.read_csv(fp_path, sep='\t')

fn['confusion_matrix'], tp['confusion_matrix'], tn['confusion_matrix'], fp['confusion_matrix'] = 'FN', 'TP', 'TN', 'FP'

tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
color_dict = snakemake.params.color_dict
color_dict['Others'] = 'lightgrey'

# Create avg_tpm and log2_avg_tpm columns in tpm_df
tpm_df['avg_tpm'] = tpm_df.filter(regex='^[A-Z].*_[123]$').mean(axis=1)
tpm_df['avg_tpm'] = tpm_df['avg_tpm'] + 0.0001  # add pseudocount to be able to compute log afterwards
tpm_df['gene_id_sno'] = tpm_df['gene_id']
tpm_df['log2_avg_tpm'] = np.log2(tpm_df['avg_tpm'])
tpm_df = tpm_df[['gene_id_sno', 'log2_avg_tpm', 'avg_tpm']]

# Merge each confusion value df to tpm_df
fn = fn.merge(tpm_df, how='left', on='gene_id_sno')
tp = tp.merge(tpm_df, how='left', on='gene_id_sno')
tn = tn.merge(tpm_df, how='left', on='gene_id_sno')
fp = fp.merge(tpm_df, how='left', on='gene_id_sno')
fn = fn[['gene_id_sno', 'log2_avg_tpm', 'avg_tpm', 'confusion_matrix']]
tp = tp[['gene_id_sno', 'log2_avg_tpm', 'avg_tpm', 'confusion_matrix']]
tn = tn[['gene_id_sno', 'log2_avg_tpm', 'avg_tpm', 'confusion_matrix']]
fp = fp[['gene_id_sno', 'log2_avg_tpm', 'avg_tpm', 'confusion_matrix']]

# Find all snoRNAs not included as a real FN, TP, TN or FP
real_confusion_val_ids = []  # all real confusion value snoRNAs (those always predicted as such across models and iterations)
for df in [fn, tp, tn, fp]:
    ids = list(df.gene_id_sno)
    real_confusion_val_ids.append(ids)
    print(len(df))
real_confusion_val_ids = [val for sublist in real_confusion_val_ids for val in sublist]  # flatten list of lists into list of values
others = tpm_df[~tpm_df.gene_id_sno.isin(real_confusion_val_ids)]
others['confusion_matrix'] = 'Others'
print(len(others))
# Create the violin plot
concat_df = pd.concat([fn, tp, tn, fp, others])

ft.violin_wo_swarm(concat_df, 'confusion_matrix', 'log2_avg_tpm', None, 'Confusion value', 'log2(average TPM)', '',
            color_dict, snakemake.output.violin)
