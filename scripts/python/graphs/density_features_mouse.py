#!/usr/bin/python3
import pandas as pd
import functions as ft
import numpy as np
df = pd.read_csv(snakemake.input.df, sep='\t')
snoRNA_type_df = pd.read_csv(snakemake.input.snoRNA_type_df, sep='\t')
sno_type = str(snakemake.wildcards.sno_type)
sno_type = sno_type[0] + '/' + sno_type[1:]

# Keep only C/D or H/ACA snoRNAs
df = df.merge(snoRNA_type_df, how='left', left_on='gene_id_sno', right_on='gene_id')
df = df[df['snoRNA_type'] == sno_type]

# Generate a density plot of numerical features with a hue of abundance_cutoff
hues = list(pd.unique(df['abundance_cutoff']))
df_list = []
colors = []
color_dict = snakemake.params.hue_color

for hue in hues:
    temp_df = df[df['abundance_cutoff'] == hue][snakemake.wildcards.mouse_numerical_features]
    df_list.append(temp_df)
    color = color_dict[hue]
    colors.append(color)




ft.density_x(df_list, snakemake.wildcards.mouse_numerical_features, 'Density', 'linear', '',
        colors, hues, snakemake.output.density_features)
