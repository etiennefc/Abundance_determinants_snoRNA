#!/usr/bin/python3
import pandas as pd
import numpy as np

""" Generate a merged dataframe of all snoRNA features and labels that will be
    used by the predictor. The labels are 'expressed' or 'not_expressed' using
    the column 'abundance_cutoff'. The features are: host_expressed, snoRNA 
    structure stability (Minimal Free Energy or MFE), snoRNA terminal stem 
    stability (MFE), snoRNA conservation and hamming_distance_box (per box or
    global hamming distance). """

# Labels (abundance_cutoff) and feature abundance_cutoff_host
tpm_df_labels = pd.read_csv(snakemake.input.abundance_cutoff, sep='\t')
tpm_df_labels = tpm_df_labels[['gene_id', 'gene_name', 'abundance_cutoff', 'abundance_cutoff_host']]
tpm_df_labels = tpm_df_labels.rename(columns={'gene_id': 'gene_id_sno'})

# Other features
sno_mfe = pd.read_csv(snakemake.input.sno_structure_mfe, sep='\t',
            names=['gene_id_sno', 'sno_mfe'])
terminal_stem_mfe = pd.read_csv(snakemake.input.terminal_stem_mfe, sep='\t',
                    names=['gene_id_sno', 'terminal_stem_mfe'])

hamming = pd.read_csv(snakemake.input.hamming_distance_box, sep='\t')
hamming = hamming[['gene_id', 'combined_box_hamming']]  # get combined hamming distance for all boxes in a snoRNA
hamming = hamming.rename(columns={'gene_id': 'gene_id_sno'})


# Merge iteratively all of these dataframes
df_list = [tpm_df_labels, hamming, sno_mfe, terminal_stem_mfe]

df_label = df_list[0]
temp = [df_label]
for i, df in enumerate(df_list[1:]):
    if i == 0:
        df_temp = temp[0].merge(df, how='left', on='gene_id_sno')
        temp.append(df_temp)
    else:
        df_temp = temp[i].merge(df, how='left', on='gene_id_sno')
        temp.append(df_temp)

final_df = temp[-1]  # get the last df in temp, i.e. the final merged df

final_df.to_csv(snakemake.output.feature_df, index=False, sep='\t')
