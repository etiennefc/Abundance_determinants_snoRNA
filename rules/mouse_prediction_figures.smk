import os


rule violin_models_accuracies_iterations_mouse:
    """ For each model (log_reg, svc and rf), represent a violin plot of the 
        accuracies across the 10 iterations based on the predictions of the 
        abundance status of mouse snoRNAs."""
    input:
        accuracies = expand(rules.test_accuracy_mouse.output.test_accuracy, models2=config['models3'], manual_iteration=config['manual_iteration'])
    output:
        violin = os.path.join(config['figures']['violin'], 'models_accuracies_iterations_mouse.svg')
    params:
        color_dict = config['colors_complex']['model_colors']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/violin_models_accuracies_iterations_mouse.py"


rule get_consensus_confusion_value_mouse:
    """ Define the consensus confusion value across the 3 models and 10 iterations. 
        Choose the confusion value based on the highest number of time predicted 
        as such across the 30 models. If 2 confusion values have an equal number of 
        votes (15 vs 15), remove randomly 1 iteration and the equality will then be 
        broken and a predominant confusion value will be chosen."""
    input:
        confusion_val_df = expand(rules.confusion_matrix_f1_mouse.output.info_df, models2=config['models3'], manual_iteration=config['manual_iteration'])
    output:
        consensus_conf_val_df = os.path.join(config['path']['confusion_matrix_f1_mouse'],
                                    'consensus_confusion_value_per_sno.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_consensus_confusion_value_mouse.py"

rule pie_confusion_values_mouse:
    """ Create a pie chart of the number of mouse snoRNAs per confusion value."""
    input:
        confusion_value_per_sno = rules.get_consensus_confusion_value_mouse.output.consensus_conf_val_df
    output:
        pie = os.path.join(config['figures']['pie'], 'sno_number_per_confusion_value_mouse.svg')
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/pie_confusion_values_mouse.py"

rule donut_label_sno_type_mouse:
    """ Generate a donut chart of the number and % of expressed vs not expressed
        snoRNAs (outer donut) and per sno_type (inner donut) for mouse snoRNAs."""
    input:
        df = rules.find_mouse_snoRNA_labels_w_length.output.tpm_label_df
    output:
        donut = os.path.join(config['figures']['donut'],
                            'abundance_status_sno_type_mouse.svg')
    params:
        label_colors = config['colors_complex']['abundance_cutoff_2'],
        sno_type_colors = config['colors_complex']['sno_type']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/donut_labels_sno_type_mouse.py"

rule donut_label_host_biotype_mouse:
    """ Generate a donut chart of the number and % of expressed vs not expressed
        snoRNAs (outer donut) and per host biotype (inner donut) for mouse snoRNAs."""
    input:
        df = rules.find_mouse_snoRNA_labels_w_length.output.tpm_label_df,
        host_biotype_df = rules.find_mouse_snoRNA_HG.output.mouse_snoRNA_HG
    output:
        donut = os.path.join(config['figures']['donut'],
                            'abundance_status_host_biotype_mouse.svg')
    params:
        label_colors = config['colors_complex']['abundance_cutoff_2'],
        host_biotype_colors = config['colors_complex']['host_biotype2']
    conda:
        "../envs/python.yaml"   
    script:
        "../scripts/python/graphs/donut_labels_host_biotype_mouse.py"


##for sno_type and host biotype
