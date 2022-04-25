import os

include: "downloads.smk"
include: "data_processing.smk"

rule add_HG_tpm_gtex:
    """ Add host gene abundance columns (for each triplicate tissue from GTEx)
        for each snoRNA within a snoRNA dataframe."""
    input:
        tpm_biotype_df = rules.add_gene_biotype.output.tpm_biotype,
        ref_HG_table = config['path']['host_gene_df'],
        gtex_tpm_df = rules.download_gtex_host_data.output.gtex_data
    output:
        sno_HG_df = config['path']['sno_tpm_df_gtex_HG']
    params:
        gtex_id_dict = config['gtex_id_tissue_dict']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/add_HG_tpm_gtex.py"

rule abundance_cutoff_gtex_HG:
    """ Define a cutoff to classify snoRNAs as 'expressed' or 'not_expressed'
        across tissues (this is the label that will be used by the predictor).
        Define the same cutoff for the host genes (used as another feature:
        HG expressed, not expressed or no HG). """
    input:
        tpm_df = rules.add_HG_tpm_gtex.output.sno_HG_df
    output:
        abundance_cutoff_df = config['path']['sno_tpm_df_cutoff_gtex_HG']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/abundance_cutoff.py"

rule add_HG_tpm_gtex_unpaired:
    """ Add host gene abundance columns (for each triplicate tissue from GTEx
        originating from unpaired tissue (not those in TGIRT-Seq datasets))
        for each snoRNA within a snoRNA dataframe."""
    input:
        tpm_biotype_df = rules.add_gene_biotype.output.tpm_biotype,
        ref_HG_table = config['path']['host_gene_df'],
        gtex_tpm_df = rules.download_gtex_host_data_unpaired.output.gtex_data
    output:
        sno_HG_df = config['path']['sno_tpm_df_gtex_HG_unpaired']
    params:
        gtex_id_dict = config['gtex_id_unpaired_tissue_dict']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/add_HG_tpm_gtex.py"

rule abundance_cutoff_gtex_HG_unpaired:
    """ Define a cutoff to classify snoRNAs as 'expressed' or 'not_expressed'
        across tissues (this is the label that will be used by the predictor).
        Define the same cutoff for the host genes (used as another feature:
        HG expressed, not expressed or no HG). """
    input:
        tpm_df = rules.add_HG_tpm_gtex_unpaired.output.sno_HG_df
    output:
        abundance_cutoff_df = config['path']['sno_tpm_df_cutoff_gtex_HG_unpaired']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/abundance_cutoff.py"
