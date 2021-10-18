#!/usr/bin/python3
import pandas as pd
import functions as ft

df = pd.read_csv(snakemake.input.all_features, sep='\t')
haca = df[df['sno_type'] == 'H/ACA']
colors_dict = snakemake.params.ab_status_color

ft.bivariate_density(haca, 'conservation_score', 'sno_mfe', 'abundance_cutoff_2',
                    snakemake.output.bivariate_density, palette=colors_dict,
                    edgecolor='grey', xlim=(-0.1,1.1), ylim=(-175,0))


expressed = haca[haca['abundance_cutoff_2'] == 'expressed']
print('Average conservation across expressed H/ACA', expressed['conservation_score'].mean())
print('Average sno_mfe across expressed H/ACA', expressed['sno_mfe'].mean())
