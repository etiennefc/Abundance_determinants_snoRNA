#!/usr/bin/python3
import pandas as pd
import functions as ft
df = pd.read_csv(snakemake.input.df, sep='\t')
feat = snakemake.wildcards.numerical_features
# Reverse the relative intron rank value (count it from the 5' end instead of from the 3' end)
if feat == 'relative_intron_rank':
    df['relative_intron_rank_switch'] = 1 - df['relative_intron_rank']
    feat = 'relative_intron_rank_switch'
print(feat)
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
    temp_df_cd = cd[cd['abundance_cutoff_2'] == hue][feat]
    df_list_cd.append(temp_df_cd)
    temp_df_haca = haca[haca['abundance_cutoff_2'] == hue][feat]
    df_list_haca.append(temp_df_haca)
    color = color_dict[hue]
    colors.append(color)
if feat == 'terminal_stem_mfe':
    ft.density_x_size(df_list_cd, feat, 'Density', 'linear', '',
            colors, hues, snakemake.output.density_features_cd, (10,8), -42, 2)
    ft.density_x_size(df_list_haca, feat, 'Density', 'linear', '',
            colors, hues, snakemake.output.density_features_haca, (10,8), -42, 2)
elif feat == 'relative_intron_rank_switch':    
    ft.density_x_size(df_list_cd, 'relative_intron_rank_switch', 'Density', 'linear', '',
            colors, hues, snakemake.output.density_features_cd, (10,8), -0.05, 1.05)
    ft.density_x_size(df_list_haca, 'relative_intron_rank_switch', 'Density', 'linear', '',
            colors, hues, snakemake.output.density_features_haca, (10,8), -0.05, 1.05)
else:
    ft.density_x(df_list_cd, feat, 'Density', 'linear', '',
            colors, hues, snakemake.output.density_features_cd)
    ft.density_x(df_list_haca, feat, 'Density', 'linear', '',
            colors, hues, snakemake.output.density_features_haca)
