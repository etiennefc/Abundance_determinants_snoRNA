#!/usr/bin/python3
import pandas as pd
import functions as ft
import math
df = pd.read_csv(snakemake.input.df, sep='\t')

# Select intronic snoRNAs and create intron subgroup (small vs long intron)
df = df[df['host_biotype2'] != 'intergenic']
df.loc[df['intron_length'] < 5000, 'intron_subgroup'] = 'small_intron'
df.loc[df['intron_length'] >= 5000, 'intron_subgroup'] = 'long_intron'

# Separate per snoRNA type
sno_type = snakemake.wildcards.sno_type
sno_type = sno_type[0] + '/' + sno_type[1:]
sno_type_df = df[df['sno_type'] == sno_type]

# Generate a density plot of various features with a hue of abundance_cutoff_2
# per snoRNA type and intron subgroup
hues = list(pd.unique(sno_type_df['abundance_cutoff_2']))
df_list_small = []
df_list_long = []
colors = []
color_dict = snakemake.params.hue_color

# Create function to get min and max values to set the xlim accordingly on the density plot x-axis
def get_min_max(df):
    if snakemake.wildcards.intron_group_feature in ['terminal_stem_mfe', 'sno_mfe']:
        min = df[snakemake.wildcards.intron_group_feature].min()
        min = -1 * (-min + (10 - -min % 10))  # get the closest negative number that is a multiple of 10 going downward in negative values
        max = 0
        min2, max2 = None, None
    elif snakemake.wildcards.intron_group_feature == 'conservation_score':
        min, max = -0.05, 1.05
        min2, max2 = None, None
    elif snakemake.wildcards.intron_group_feature == 'dist_to_bp':  # these values vary too much between intron subgroups to have the same xlims
        min = 0                                                     # this is why we return 2 min and 2 max values (one for each intron subgroup)
        df_small, df_long = df[df['intron_subgroup'] == 'small_intron'], df[df['intron_subgroup'] == 'long_intron']
        max = df_small[snakemake.wildcards.intron_group_feature].max()
        max = max + (10 - max % 10)  # get the closest positive number that is a multiple of 10 going upward in positive values
        min2 = 0
        max2 = df_long[snakemake.wildcards.intron_group_feature].max()
        max2 = max2 + (10 - max2 % 10)  # get the closest positive number that is a multiple of 10 going upward in positive values
    else:
        min = 0
        max = df[snakemake.wildcards.intron_group_feature].max()
        max = max + (10 - max % 10)  # get the closest positive number that is a multiple of 10 going upward in positive values
        min2, max2 = None, None

    return min, max, min2, max2

# Separate per intron subgroup and get xlims for the density plots
small_df = sno_type_df[sno_type_df['intron_subgroup'] == 'small_intron']
long_df = sno_type_df[sno_type_df['intron_subgroup'] == 'long_intron']
min, max, dist_to_bp_min2, dist_to_bp_max2 = get_min_max(sno_type_df)

for hue in hues:
    temp_df_small = small_df[small_df['abundance_cutoff_2'] == hue][snakemake.wildcards.intron_group_feature]
    temp_df_long = long_df[long_df['abundance_cutoff_2'] == hue][snakemake.wildcards.intron_group_feature]
    df_list_small.append(temp_df_small)
    df_list_long.append(temp_df_long)
    colors.append(color_dict[hue])

if snakemake.wildcards.intron_group_feature == 'dist_to_bp':
    ft.density_x_size(df_list_small, snakemake.wildcards.intron_group_feature, 'Density', 'linear', '',
            colors, hues, snakemake.output.density_small, (12, 6), min, max)
    ft.density_x_size(df_list_long, snakemake.wildcards.intron_group_feature, 'Density', 'linear', '',
            colors, hues, snakemake.output.density_long, (12, 6), dist_to_bp_min2, dist_to_bp_max2)
else:
    if (math.isnan(max) == True) | (math.isnan(min) == True):  # when looking at C,D,C',D' hamming for H/ACA snoRNAs or H,ACA hamming for C/D, the result is NaN values
        min_modified, max_modified = 0, 1
        ft.density_x_size(df_list_small, snakemake.wildcards.intron_group_feature, 'Density', 'linear', '',
                colors, hues, snakemake.output.density_small, (12, 6), min_modified, max_modified)
        ft.density_x_size(df_list_long, snakemake.wildcards.intron_group_feature, 'Density', 'linear', '',
                colors, hues, snakemake.output.density_long, (12, 6), min_modified, max_modified)
    else:
        ft.density_x_size(df_list_small, snakemake.wildcards.intron_group_feature, 'Density', 'linear', '',
                colors, hues, snakemake.output.density_small, (12, 6), min, max)
        ft.density_x_size(df_list_long, snakemake.wildcards.intron_group_feature, 'Density', 'linear', '',
                colors, hues, snakemake.output.density_long, (12, 6), min, max)
