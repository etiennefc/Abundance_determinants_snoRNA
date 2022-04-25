#!/usr/bin/python3
import pandas as pd
import functions as ft

df = pd.read_csv(snakemake.input.df, sep='\t')
snoRNA_type_df = pd.read_csv(snakemake.input.snoRNA_type_df, sep='\t')
sno_type = str(snakemake.wildcards.sno_type)
sno_type = sno_type[0] + '/' + sno_type[1:]

# Keep only C/D or H/ACA snoRNAs
df = df.merge(snoRNA_type_df, how='left', left_on='gene_id_sno', right_on='gene_id')
df = df[df['snoRNA_type'] == sno_type]

# Generate a bar chart of categorical features with a hue of abundance_cutoff_2
# per sno_type
counts_per_feature = ft.count_list_x(df, 'abundance_cutoff',
                    list(snakemake.params.hue_color.keys()),
                    'abundance_cutoff_host')
percent = ft.percent_count(counts_per_feature)

ft.stacked_bar(percent, sorted(list(df['abundance_cutoff'].unique())),
                list(snakemake.params.hue_color.keys()), '', 'Abundance status of snoRNAs',
                'Proportion of snoRNAs (%)', snakemake.params.hue_color,
            snakemake.output.bar)
