#!/usr/bin/python3
import pandas as pd

""" Keep only one of the top predicting feature from the orginal
    dataset. Either conservation_score_norm,
    terminal_stem_mfe_norm, sno_mfe_norm, host_expressed, intergenic or
    dist_to_bp_norm."""

df = pd.read_csv(snakemake.input.df, sep='\t')

if snakemake.wildcards.one_feature == "all_four":
    df = df[['gene_id_sno', "conservation_score_norm", "terminal_stem_mfe_norm",
            "sno_mfe_norm", "host_expressed", 'label']]
elif snakemake.wildcards.one_feature == "all_five":
    df = df[['gene_id_sno', "conservation_score_norm", "terminal_stem_mfe_norm",
            "sno_mfe_norm", "host_expressed", "intron_number_norm", 'label']]
else:
    df = df[['gene_id_sno', snakemake.wildcards.one_feature, 'label']]

df.to_csv(snakemake.output.df_one_feature, index=False, sep='\t')
