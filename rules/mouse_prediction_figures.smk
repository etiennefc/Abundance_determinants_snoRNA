import os

include: "cv_train_test_for_species_prediction_top4.smk"
include: "cv_train_test_for_species_prediction_top4_random_state.smk"
include: "cv_train_test_for_species_prediction_top4_log_reg_thresh.smk"

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

rule get_consensus_confusion_value_per_model_mouse:
    """ Define the consensus confusion value across the 10 iterations for each model.
        Choose the confusion value based on the highest number of time predicted
        as such across the 10 models. If 2 confusion values have an equal number of
        votes (5 vs 5), remove randomly 1 iteration and the equality will then be
        broken and a predominant confusion value will be chosen."""
    input:
        confusion_val_df = expand(rules.confusion_matrix_f1_mouse.output.info_df, manual_iteration=config['manual_iteration'], allow_missing=True)
    output:
        consensus_conf_val_df = os.path.join(config['path']['confusion_matrix_f1_mouse'],
                                    'consensus_confusion_value_per_sno_{models2}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_consensus_confusion_value_per_model_mouse.py"

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

rule pie_confusion_values_per_model_mouse:
    """ Create a pie chart of the number of mouse snoRNAs per confusion value
        per model (svc, log_reg and rf)."""
    input:
        confusion_value_per_sno = rules.get_consensus_confusion_value_per_model_mouse.output.consensus_conf_val_df
    output:
        pie = os.path.join(config['figures']['pie'], 'sno_number_per_confusion_value_mouse_{models2}.svg')
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/pie_confusion_values_mouse.py"

rule pie_confusion_values_per_model_species_prediction:
    """ Create a pie chart of the number of mouse snoRNAs per confusion value
        per model (svc, log_reg and rf) trained with all human snoRNAs."""
    input:
        confusion_value_per_sno = rules.confusion_matrix_f1_species_prediction_top4_rs.output.info_df
    output:
        pie = os.path.join(config['figures']['pie'], 'sno_number_per_confusion_value_species_prediction_{models2}_{rs}.svg')
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/pie_confusion_values_species_prediction.py"

rule pie_confusion_values_species_prediction_log_reg_thresh:
    """ Create a pie chart of the number of mouse snoRNAs per confusion value
        for the log_reg thresh model trained with all human snoRNAs."""
    input:
        confusion_value_per_sno = expand(rules.confusion_matrix_f1_species_prediction_top4_log_reg_thresh.output.info_df, rs="42")
    output:
        pie = os.path.join(config['figures']['pie'], 'sno_number_per_confusion_value_species_prediction_log_reg_thresh.svg')
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/pie_confusion_values_species_prediction.py"

rule donut_confusion_values_host_biotype_species_prediction_log_reg_thresh:
    """ Generate a donut chart of the number and % of confusion value
        snoRNAs (outer donut) and per host biotype (inner donut) for mouse
        snoRNAs using the log_reg_thresh model."""
    input:
        host_biotype_df = rules.find_mouse_snoRNA_HG.output.mouse_snoRNA_HG,
        confusion_value_per_sno = expand(rules.confusion_matrix_f1_species_prediction_top4_log_reg_thresh.output.info_df, rs="42")
    output:
        donut = os.path.join(config['figures']['donut'],
                            'confusion_value_host_biotype_mouse.svg')
    params:
        conf_val_colors = config['colors_complex']['confusion_value'],
        host_biotype_colors = config['colors_complex']['host_biotype2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/donut_confusion_values_host_biotype_species_prediction_log_reg_thresh.py"


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

rule violin_tpm_confusion_value_mouse:
    """ Create a violin plot of the log2(TPM) for all confusion value snoRNAs."""
    input:
        sno_per_confusion_value = rules.get_consensus_confusion_value_mouse.output.consensus_conf_val_df,
        tpm_df = rules.find_mouse_snoRNA_labels_w_length.output.tpm_label_df
    output:
        violin = os.path.join(config['figures']['violin'], "avg_tpm_per_confusion_value_mouse.svg")
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/violin_tpm_confusion_value_mouse.py"

rule violin_tpm_confusion_value_per_model_mouse:
    """ Create a violin plot of the log2(TPM) for all confusion value snoRNAs per model (log_reg, svc and rf)."""
    input:
        sno_per_confusion_value = rules.get_consensus_confusion_value_per_model_mouse.output.consensus_conf_val_df,
        tpm_df = rules.find_mouse_snoRNA_labels_w_length.output.tpm_label_df
    output:
        violin = os.path.join(config['figures']['violin'], "avg_tpm_per_confusion_value_mouse_{models2}.svg")
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/violin_tpm_confusion_value_mouse.py"

rule violin_tpm_confusion_value_log_reg_thresh_mouse:
    """ Create a violin plot of the log2(TPM) for all confusion value snoRNAs for the log_reg thresh model."""
    input:
        sno_per_confusion_value = expand(rules.confusion_matrix_f1_species_prediction_top4_log_reg_thresh.output.info_df, rs="42"),
        tpm_df = rules.find_mouse_snoRNA_labels_w_length.output.tpm_label_df
    output:
        violin = os.path.join(config['figures']['violin'], "avg_tpm_per_confusion_value_mouse_log_reg_thresh.svg")
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/violin_tpm_confusion_value_log_reg_thresh_mouse.py"

rule scatter_accuracies_species_prediction_top4:
    """ Generate a connected scatter plot for each model to show their accuracy
        of prediction on the CV, training and test sets to highlight possible
        overfitting."""
        input:
            cv_accuracy = expand(os.path.join(config['path']['hyperparameter_tuning'],
                                    '{models2}_best_params_top4_species_prediction.tsv'),
                                    **config),
            training_accuracy = expand(os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_top4_species_prediction.tsv'),
                                    **config),
            test_accuracy = expand(os.path.join(config['path']['test_accuracy_mouse'],
                                    '{models2}_test_accuracy_top4_species_prediction.tsv'),
                                    **config),
        output:
            scatter = os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_species_prediction.svg')
        params:
            colors = config['colors_complex']['model_colors']
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/python/graphs/scatter_accuracies.py"

rule scatter_accuracies_species_prediction_top4_random_state:
    """ Generate a connected scatter plot for each model to show their accuracy
        of prediction on the CV, training and test sets to highlight possible
        overfitting."""
        input:
            cv_accuracy = expand(os.path.join(config['path']['hyperparameter_tuning'],
                                    '{models2}_best_params_top4_species_prediction_{rs}.tsv'),
                                    **config),
            training_accuracy = expand(os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_top4_species_prediction_{rs}.tsv'),
                                    **config),
            test_accuracy = expand(os.path.join(config['path']['test_accuracy_mouse'],
                                    '{models2}_test_accuracy_top4_species_prediction_{rs}.tsv'),
                                    **config),
        output:
            scatter = os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_species_prediction_rs.svg')
        params:
            colors = config['colors_complex']['model_colors']
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/python/graphs/scatter_accuracies_10_iterations.py"

rule scatter_accuracies_species_prediction_top4_random_state_w_log_reg_thresh:
    """ Generate a connected scatter plot for each model to show their accuracy
        of prediction on the CV, training and test sets to highlight possible
        overfitting."""
        input:
            cv_accuracy = expand(os.path.join(config['path']['hyperparameter_tuning'],
                                    '{models2}_best_params_top4_species_prediction_{rs}.tsv'),
                                    **config),
            training_accuracy = expand(os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_top4_species_prediction_{rs}.tsv'),
                                    models2=["svc", "rf", "gbm", "knn"], rs=config['rs']),
            training_accuracy_log_reg_thresh = expand(rules.train_test_accuracy_species_prediction_top4_log_reg_thresh.output.training_accuracy, **config),
            test_accuracy = expand(os.path.join(config['path']['test_accuracy_mouse'],
                                    '{models2}_test_accuracy_top4_species_prediction_{rs}.tsv'),
                                    models2=["svc", "rf", "gbm", "knn"], rs=config['rs']),
            test_accuracy_log_reg_thresh = expand(rules.train_test_accuracy_species_prediction_top4_log_reg_thresh.output.test_accuracy, **config)
        output:
            scatter = os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_species_prediction_rs_w_log_reg_thresh.svg')
        params:
            colors = config['colors_complex']['model_colors']
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/python/graphs/scatter_accuracies_species_prediction_top4_random_state_w_log_reg_thresh.py"

rule scatter_accuracies_species_prediction_top4_w_log_reg_thresh:
    """ Generate a connected scatter plot for each model to show their accuracy
        of prediction on the CV, training and test sets to highlight possible
        overfitting (including the log_reg thresh model)."""
        input:
            cv_accuracy = expand(os.path.join(config['path']['hyperparameter_tuning'],
                                    '{models2}_best_params_top4_species_prediction.tsv'),
                                    **config),
            training_accuracy = expand(os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_top4_species_prediction.tsv'),
                                    models2=["svc", "rf", "gbm", "knn"]),
            training_accuracy_log_reg_thresh = expand(rules.train_test_accuracy_species_prediction_top4_log_reg_thresh.output.training_accuracy, rs="42"),
            test_accuracy = expand(os.path.join(config['path']['test_accuracy_mouse'],
                                    '{models2}_test_accuracy_top4_species_prediction.tsv'),
                                    models2=["svc", "rf", "gbm", "knn"]),
            test_accuracy_log_reg_thresh = expand(rules.train_test_accuracy_species_prediction_top4_log_reg_thresh.output.test_accuracy, rs="42")
        output:
            scatter = os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_species_prediction_w_log_reg_thresh.svg')
        params:
            colors = config['colors_complex']['model_colors']
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/python/graphs/scatter_accuracies_species_prediction_top4_w_log_reg_thresh.py"
