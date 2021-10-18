#!/usr/bin/python3
import pandas as pd
import functions as ft

df = pd.read_csv(snakemake.input.rank_features_df, sep='\t')
all_features = pd.read_csv(snakemake.input.all_features_df, sep='\t')

# Get top 5 features across models
top_5 = df.groupby('feature')['feature_rank'].median().sort_values(ascending=True).index.to_list()[0:5]
top_5 = ' '.join(top_5).replace('host_expressed_norm', 'host_expressed').split()
top_5.extend(['C/D', 'label'])
top_5_df = all_features.filter(top_5, axis=1)

# Create df for C/D and H/ACA snoRNAs
cd = top_5_df[top_5_df['C/D'] == 1]
cd = cd.drop('C/D', axis=1)
haca = top_5_df[top_5_df['C/D'] == 0]
haca = haca.drop('C/D', axis=1)

# Replace string labels by numeric labels
hue_color = snakemake.params.hue_color
hue_color[0] = hue_color.pop('not_expressed')
hue_color[1] = hue_color.pop('expressed')

# Generate a pairplot of top 5 features across models with a hue of label
# for either C/D and H/ACA snoRNAs separately
ft.pairplot(cd, 'label', snakemake.params.hue_color,
            snakemake.output.pairplot_cd)
ft.pairplot(haca, 'label', snakemake.params.hue_color,
        snakemake.output.pairplot_haca)
