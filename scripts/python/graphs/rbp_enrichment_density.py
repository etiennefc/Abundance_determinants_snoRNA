#!/usr/bin/python3
import pandas as pd
import functions as ft
import numpy as np
from functools import reduce

feature_df = pd.read_csv(snakemake.input.all_feature_df, sep='\t')
rbp_bed_paths = snakemake.input.beds
output_simple = snakemake.output.density
color_dict = snakemake.params.color_dict
colors = [color_dict['expressed'], color_dict['not_expressed']]
bed_dfs = []
for i, path in enumerate(rbp_bed_paths):
    rbp = path.split('/')[-1].split('_mapped')[0]
    print(rbp)
    bed_df = pd.read_csv(path, sep='\t', names=["chr", "start", "end", "gene_id_sno", "dot", "strand", f"{rbp}_score"])
    bed_df = bed_df[['gene_id_sno', f'{rbp}_score']]
    bed_df[f'{rbp}_score'] = bed_df[f'{rbp}_score'].replace('.', 0).astype(float)
    bed_df[f'{rbp}_score'] = bed_df[f'{rbp}_score'] + 0.000001  # add pseudocount so that a log can be computed
    bed_df = bed_df.merge(feature_df, how='left', on='gene_id_sno')
    bed_df[f'{rbp}_score_log10'] = np.log10(bed_df[f'{rbp}_score'])

    # Get only intergenic or intronic snoRNAs
    #bed_df = bed_df[bed_df['abundance_cutoff_host'] == 'intergenic']
    # Create a density plot for each RBP enrichment to compare expressed vs not expressed snoRNAs
    expressed, not_expressed = bed_df[bed_df['abundance_cutoff_2'] == 'expressed'], bed_df[bed_df['abundance_cutoff_2'] == 'not_expressed']
    ft.density_x([expressed[f'{rbp}_score_log10'], not_expressed[f'{rbp}_score_log10']], f'{rbp} enrichment score (log10)', 'Density', 'linear', '', colors,
                ['expressed', 'not_expressed'], output_simple[i])
    bed_df = bed_df.drop(columns='abundance_cutoff_2')
    bed_dfs.append(bed_df)


# Create a combined RBP enrichment score
df_merged = reduce(lambda left,right: pd.merge(left,right,on=['gene_id_sno'],
                                            how='left'), bed_dfs)
df_merged = df_merged.filter(regex='(gene_id_sno|_score$)')
df_merged = df_merged.drop(columns='FBL_mnase_score')
df_merged['combined_rbp_score'] = df_merged.filter(regex='_score$').sum(axis=1)
df_merged = df_merged.merge(feature_df[['gene_id_sno', 'abundance_cutoff_2']], on='gene_id_sno', how='left')
df_merged['combined_rbp_score_log10'] = np.log10(df_merged['combined_rbp_score'])
df_merged.to_csv(snakemake.output.combined_rbp_score_df, sep='\t', index=False)


expressed_merged, not_expressed_merged = df_merged[df_merged['abundance_cutoff_2'] == 'expressed'], df_merged[df_merged['abundance_cutoff_2'] == 'not_expressed']

# Create a density plot for the combined RBP enrichment to compare expressed vs not expressed snoRNAs
ft.density_x([expressed_merged['combined_rbp_score_log10'], not_expressed_merged['combined_rbp_score_log10']], 'Combined RBP enrichment score (log10)', 'Density', 'linear', '', colors,
            ['expressed', 'not_expressed'], snakemake.output.density_combined)
