#!/usr/bin/python3
import pandas as pd
import functions as ft
df = pd.read_csv(snakemake.input.df, sep='\t')
cd = df[df['sno_type'] == 'C/D']
haca = df[df['sno_type'] == 'H/ACA']

# Generate a density plot of numerical features with a hue of abundance_cutoff_2
# for either C/D and H/ACA snoRNAs separately
hues = list(pd.unique(df['abundance_cutoff_2']))
df_list_cd = []
df_list_haca = []
colors = []
color_dict = snakemake.params.hue_color

for hue in hues:
    temp_df_cd = cd[cd['abundance_cutoff_2'] == hue][snakemake.wildcards.numerical_features_scaled]
    df_list_cd.append(temp_df_cd)
    temp_df_haca = haca[haca['abundance_cutoff_2'] == hue][snakemake.wildcards.numerical_features_scaled]
    df_list_haca.append(temp_df_haca)
    color = color_dict[hue]
    colors.append(color)


ft.density_x(df_list_cd, snakemake.wildcards.numerical_features_scaled, 'Density', 'linear', '',
        colors, hues, snakemake.output.density_features_cd)

ft.density_x(df_list_haca, snakemake.wildcards.numerical_features_scaled, 'Density', 'linear', '',
        colors, hues, snakemake.output.density_features_haca)
