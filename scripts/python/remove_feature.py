#!/usr/bin/python3
import pandas as pd

""" Remove the top predicting feature from the orginal
    dataset containing all snoRNA features. Either conservation_score_norm,
    terminal_stem_mfe_norm, sno_mfe_norm, host_expressed or all_four features."""

df = pd.read_csv(snakemake.input.df, sep='\t')

if snakemake.wildcards.feature_effect == "all_four":
    df = df.drop(columns=["conservation_score_norm", "terminal_stem_mfe_norm",
                        "sno_mfe_norm", "host_expressed"])
elif snakemake.wildcards.feature_effect == "top_10":
    df = df.drop(columns=["conservation_score_norm", "terminal_stem_mfe_norm",
                        "sno_mfe_norm", "host_expressed", "intron_number_norm",
                        "dist_to_bp_norm", "intron_length_norm", "sno_length_norm", "intergenic.3",
                        "rRNA", "Orphan", "distance_downstream_exon_norm"])
elif snakemake.wildcards.feature_effect == "top_10_intergenic":
    df = df.drop(columns=["conservation_score_norm", "terminal_stem_mfe_norm",
                        "sno_mfe_norm", "host_expressed", "intron_number_norm",
                        "dist_to_bp_norm", "intron_length_norm", "sno_length_norm", "intergenic.3",
                        "rRNA", "Orphan", "distance_downstream_exon_norm", "intergenic", "intergenic.1", "intergenic.2", "intergenic.4"])
elif snakemake.wildcards.feature_effect == "top_10_all":
    df = df.drop(columns=["conservation_score_norm", "terminal_stem_mfe_norm",
                        "sno_mfe_norm", "host_expressed", "intron_number_norm",
                        "dist_to_bp_norm", "intron_length_norm", "sno_length_norm", "intergenic.3",
                        "rRNA", "Orphan", "distance_downstream_exon_norm", "intergenic", "intergenic.1",
                        "intergenic.2", "intergenic.4", "host_not_expressed", "snRNA", "non_coding",
                        "protein_coding","False", "True", "dual_initiation", "simple_initiation"])
elif snakemake.wildcards.feature_effect == "top_11_all":
    df = df.drop(columns=["conservation_score_norm", "terminal_stem_mfe_norm",
                        "sno_mfe_norm", "host_expressed", "intron_number_norm",
                        "dist_to_bp_norm", "intron_length_norm", "sno_length_norm", "intergenic.3",
                        "rRNA", "Orphan", "distance_downstream_exon_norm", "intergenic", "intergenic.1",
                        "intergenic.2", "intergenic.4", "host_not_expressed", "snRNA", "non_coding",
                        "protein_coding","False", "True", "dual_initiation", "simple_initiation", "distance_upstream_exon_norm"])                        
elif snakemake.wildcards.feature_effect == "Other":
    df = df[["gene_id_sno", "Other", "label"]]


### Numerical vs categorical feature
###intrinsinc vs extrinsinc features


else:
    df = df.drop(columns=[snakemake.wildcards.feature_effect])

df.to_csv(snakemake.output.df_wo_feature, index=False, sep='\t')
