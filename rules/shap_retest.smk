import os

## Create SHAP explanations on interesting data points not in the test set, but in the training set

rule get_variable_sno_label_same_host:
    """ Extract snoRNAs that are encoded within the same host gene but that
        differ in their label (at least one snoRNA that is expressed and one
        that is not expressed within the same host gene)."""
    input:
        tpm_df = config['path']['sno_tpm_df_cutoff']
    output:
        variable_sno_labels_df = config['path']['variable_sno_label_same_host']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_variable_sno_label_same_host.py"
