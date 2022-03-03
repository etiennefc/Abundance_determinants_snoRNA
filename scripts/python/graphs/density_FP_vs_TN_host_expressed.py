#!/usr/bin/python3
import pandas as pd
import functions as ft
import numpy as np

feature = snakemake.wildcards.top_10_numerical_features
color_dict = snakemake.params.color_dict
colors = [color_dict['FP'], color_dict['TN']]
sno_per_confusion_value = snakemake.input.sno_per_confusion_value
feature_df = pd.read_csv(snakemake.input.all_features_df, sep='\t')
fp = pd.read_csv([path for path in sno_per_confusion_value if 'FP' in path][0], sep='\t')
tn = pd.read_csv([path for path in sno_per_confusion_value if 'TN' in path][0], sep='\t')

# Select only real confusion value (those always predicted as such across iterations and models)
real_fp = list(set(fp.gene_id_sno.to_list()) - set(tn.gene_id_sno.to_list()))
real_tn = list(set(tn.gene_id_sno.to_list()) - set(fp.gene_id_sno.to_list()))
real_fp_df = feature_df[feature_df['gene_id_sno'].isin(real_fp)].drop_duplicates()
real_tn_df = feature_df[feature_df['gene_id_sno'].isin(real_tn)].drop_duplicates()

# Select only TN that are within an expressed HG
real_tn_df = real_tn_df[real_tn_df['abundance_cutoff_host'] == 'host_expressed']

ft.density_x([real_fp_df[feature], real_tn_df[feature]], feature, 'Density', 'linear',
            'FP vs TN that have an expressed HG', colors, ['FP', 'TN with an expressed HG'], snakemake.output.density)
