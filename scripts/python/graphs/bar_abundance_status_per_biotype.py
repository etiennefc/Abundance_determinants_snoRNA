#!/usr/bin/python3
import pandas as pd
import functions as ft

df = pd.read_csv(snakemake.input.abundance_cutoff_df, sep='\t')
print(df)
biotypes = ['protein_coding', 'snRNA', 'snoRNA', 'tRNA', 'lncRNA']

# Keep only protein_coding, lncRNA, snRNA, snoRNA and tRNA
df = df[df['gene_biotype'].isin(biotypes)]
print(df)


# Generate a bar chart of categorical features with a hue of gene_biotype
counts_per_feature = ft.count_list_x(df, 'gene_biotype',
                    list(snakemake.params.colors.keys()),
                    'abundance_cutoff_2')
percent = ft.percent_count(counts_per_feature)

ft.stacked_bar(percent, sorted(list(df['gene_biotype'].unique())),
                list(snakemake.params.colors.keys()), '', 'Biotype',
                'Proportion of RNAs (%)', snakemake.params.colors,
            snakemake.output.bar)
