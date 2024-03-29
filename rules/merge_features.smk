import os

include: "gtex_HG_cutoff.smk"

rule merge_feature_df:
    """ Merge all feature columns (numerical or categorical) inside one
        dataframe."""
    input:
        abundance_cutoff = rules.abundance_cutoff.output.abundance_cutoff_df,
        sno_length = rules.sno_length.output.sno_length,
        snodb_nmd_di_promoters = rules.format_snodb.output.snodb_formatted,
        location_and_branchpoint = rules.get_best_bp.output.sno_distance_bp,
        sno_structure_mfe = rules.fasta_to_tsv.output.mfe_final,
        terminal_stem_mfe = rules.fasta_to_tsv_terminal_stem_mfe.output.mfe_stem_final,
        terminal_stem_length_score = rules.get_terminal_stem_length.output.length_stem,
        hamming_distance_box = rules.hamming_distance_box_all.output.hamming_distance_box_df
    output:
        feature_df = config['path']['feature_df']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/merge_features.py"

rule merge_feature_df_gtex_HG:
    """ Merge all feature columns (numerical or categorical) inside one
        dataframe. Use the host_abundance_cutoff from GTEx data."""
    input:
        abundance_cutoff = rules.abundance_cutoff_gtex_HG.output.abundance_cutoff_df,
        sno_length = rules.sno_length.output.sno_length,
        snodb_nmd_di_promoters = rules.format_snodb.output.snodb_formatted,
        location_and_branchpoint = rules.get_best_bp.output.sno_distance_bp,
        sno_structure_mfe = rules.fasta_to_tsv.output.mfe_final,
        terminal_stem_mfe = rules.fasta_to_tsv_terminal_stem_mfe.output.mfe_stem_final,
        terminal_stem_length_score = rules.get_terminal_stem_length.output.length_stem,
        hamming_distance_box = rules.hamming_distance_box_all.output.hamming_distance_box_df
    output:
        feature_df = config['path']['feature_df_gtex_HG']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/merge_features.py"

rule merge_feature_df_gtex_HG_unpaired:
    """ Merge all feature columns (numerical or categorical) inside one
        dataframe. Use the host_abundance_cutoff from GTEx data unpaired tissues."""
    input:
        abundance_cutoff = rules.abundance_cutoff_gtex_HG_unpaired.output.abundance_cutoff_df,
        sno_length = rules.sno_length.output.sno_length,
        snodb_nmd_di_promoters = rules.format_snodb.output.snodb_formatted,
        location_and_branchpoint = rules.get_best_bp.output.sno_distance_bp,
        sno_structure_mfe = rules.fasta_to_tsv.output.mfe_final,
        terminal_stem_mfe = rules.fasta_to_tsv_terminal_stem_mfe.output.mfe_stem_final,
        terminal_stem_length_score = rules.get_terminal_stem_length.output.length_stem,
        hamming_distance_box = rules.hamming_distance_box_all.output.hamming_distance_box_df
    output:
        feature_df = config['path']['feature_df_gtex_HG_unpaired']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/merge_features.py"
