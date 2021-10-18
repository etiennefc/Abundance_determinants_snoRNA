#!/usr/bin/python3
import pandas as pd
import functions as ft

df = pd.read_csv(snakemake.input.df, sep='\t')
cd = df[df['sno_type'] == 'C/D']
haca = df[df['sno_type'] == 'H/ACA']

# Generate a pairplot of numerical features with a hue of abundance_cutoff_2
# for all snoRNAs
ft.pairplot(df, 'abundance_cutoff_2', snakemake.params.hue_color,
            snakemake.output.pairplot)

# Generate a pairplot of numerical features with a hue of abundance_cutoff_2
# for either C/D and H/ACA snoRNAs separately
ft.pairplot(cd, 'abundance_cutoff_2', snakemake.params.hue_color,
            snakemake.output.pairplot_cd)
ft.pairplot(haca, 'abundance_cutoff_2', snakemake.params.hue_color,
        snakemake.output.pairplot_haca)        
