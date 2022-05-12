#!/usr/bin/python3
import pandas as pd
import numpy as np

""" Generate a merged dataframe of all snoRNA features that will be
    used by the predictor. The features are: host_expressed, snoRNA
    structure stability (Minimal Free Energy or MFE), snoRNA terminal stem
    stability (MFE), snoRNA combined box hamming. """

# abundance_cutoff_host df
host_cutoff_df = pd.read_csv(snakemake.input.abundance_cutoff, sep='\t')
host_cutoff_df = host_cutoff_df[['gene_id', 'gene_name','abundance_cutoff_host']]
host_cutoff_df = host_cutoff_df.rename(columns={'gene_id': 'gene_id_sno'})

# Other features
sno_mfe = pd.read_csv(snakemake.input.sno_structure_mfe, sep='\t',
            names=['gene_id_sno', 'sno_mfe'])
terminal_stem_mfe = pd.read_csv(snakemake.input.terminal_stem_mfe, sep='\t',
                    names=['gene_id_sno', 'terminal_stem_mfe'])

hamming = pd.read_csv(snakemake.input.hamming_distance_box, sep='\t')
hamming = hamming[['gene_id', 'combined_box_hamming']]  # get combined hamming distance for all boxes in a snoRNA
hamming = hamming.rename(columns={'gene_id': 'gene_id_sno'})


# Merge iteratively all of these dataframes
df_list = [host_cutoff_df, hamming, sno_mfe, terminal_stem_mfe]

temp = [df_list[0]]
for i, df in enumerate(df_list[1:]):
    if i == 0:
        df_temp = temp[0].merge(df, how='left', on='gene_id_sno')
        temp.append(df_temp)
    else:
        df_temp = temp[i].merge(df, how='left', on='gene_id_sno')
        temp.append(df_temp)

final_df = temp[-1]  # get the last df in temp, i.e. the final merged df

# Fill NA terminal stem value with 0 (if a stem couldn't be computed because the snoRNA is at the start/end of a chr)
final_df['terminal_stem_mfe'] = final_df['terminal_stem_mfe'].fillna(0)

final_df.to_csv(snakemake.output.feature_df, index=False, sep='\t')
