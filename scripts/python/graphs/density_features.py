#!/usr/bin/python3
import pandas as pd
import functions as ft
import numpy as np
df = pd.read_csv(snakemake.input.df, sep='\t')


# Generate a density plot of numerical features with a hue of sno_type,
# abundance_cutoff or abundance_cutoff_2
hues = list(pd.unique(df[snakemake.wildcards.feature_hue]))
df_list = []
colors = []
color_dict = snakemake.params.hue_color

#logscale_features = ['distance_upstream_exon', 'distance_downstream_exon',
#                    'dist_to_bp', 'intron_length', 'sno_length', 'sno_mfe']
logscale_features = []

print(snakemake.wildcards.numerical_features)
for hue in hues:
    if snakemake.wildcards.numerical_features in logscale_features:
        temp_df = df[df[snakemake.wildcards.feature_hue] == hue][snakemake.wildcards.numerical_features]
        temp_df['temp'] = np.log10(temp_df[snakemake.wildcards.numerical_features])
        temp_df = temp_df.drop([snakemake.wildcards.numerical_features], axis=1)
        temp_df.columns = snakemake.wildcards.numerical_features
        print(temp_df)
        df_list.append(temp_df)
    else:
        temp_df = df[df[snakemake.wildcards.feature_hue] == hue][snakemake.wildcards.numerical_features]
        df_list.append(temp_df)
    color = color_dict[hue]
    colors.append(color)




ft.density_x(df_list, snakemake.wildcards.numerical_features, 'Density', 'linear', '',
        colors, hues, snakemake.output.density_features)
