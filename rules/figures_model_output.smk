import os

include: "cv_train_test.smk"
include: "cv_train_test_wo_clusters.smk"
include: "cv_train_test_wo_feature.smk"
include: "cv_train_test_one_feature.smk"
include: "cv_train_test_scale_after_split.smk"
include: "terminal_stem.smk"
include: "structure.smk"
include: "cv_train_test_10_iterations.smk"
include: "cv_train_test_10_iterations_only_hamming.smk"
include: "cv_train_test_manual_split.smk"

rule modify_shap:
    """ Modify SHAP summary plot script (within _beewswarm.py) so that it sorts
        features vertically according to the median of abs(shap_values) and not
        by sum of shap_values."""
    output:
        fake_log = "log/modify_shap.log"
    params:
        path = config['path']['shap_path']
    shell:
        "mkdir -p log/ &&"
        "sed -i -E 's/np.sum\(np.abs/np.median\(np.abs/g' {params.path} && "
        "sed -i -E 's/global_shap_values = np.abs\(shap_values\).mean\(0\)/global_shap_values = np.median\(np.abs\(shap_values\), axis=0\)/g' {params.path}"
        "&> {output.fake_log}"

rule create_local_env_2:
    """ Import matplotlib, seaborn, sklearn, shap in local snakemake
        environment. Note that the log will be empty, it serves only to link
        this rule to a decoy output. You must answer yes to the following 'Proceed'
        statements."""
    output:
        log_create_env = "log/create_local_env_2.log"
    shell:
        """
        conda install -c conda-forge --force-reinstall seaborn
        conda install -c conda-forge --force-reinstall matplotlib
        conda install -c conda-forge --force-reinstall scikit-learn
        conda install -c conda-forge --force-reinstall shap
        conda install -c conda-forge --force-reinstall upsetplot
        conda install -c bioconda --force-reinstall viennarna
        conda install -c conda-forge --force-reinstall scipy
        &> {output}
        """

rule scatter_accuracies:
    """ Generate a connected scatter plot for each model to show their accuracy
        of prediction on the CV, training and test sets to highlight possible
        overfitting."""
        input:
            cv_accuracy = expand(os.path.join(config['path']['hyperparameter_tuning'],
                                    '{mod}_best_params_scale_after_split.tsv'),
                                    mod=['log_reg', 'svc', 'rf', 'gbm', 'knn']),
            training_accuracy = expand(os.path.join(config['path']['training_accuracy'],
                                    '{mod}_training_accuracy_scale_after_split.tsv'),
                                    mod=['log_reg', 'svc', 'rf', 'gbm', 'knn']),
            test_accuracy = expand(os.path.join(config['path']['test_accuracy'],
                                    '{mod}_test_accuracy_scale_after_split.tsv'), mod=['log_reg', 'svc', 'rf', 'gbm', 'knn']),
        output:
            scatter = os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test.svg')
        params:
            colors = config['colors_complex']['model_colors']
        script:
            "../scripts/python/graphs/scatter_accuracies.py"

rule scatter_accuracies_10_iterations:
    """ Generate a connected scatter plot for each model to show their accuracy
        of prediction on the CV, training and test sets to highlight possible
        overfitting. Show the average accuracy plus std deviation across all 10
        iterations."""
        input:
            cv_accuracy = expand(rules.hyperparameter_tuning_cv_scale_after_split_10_iterations.output.best_hyperparameters, **config),
            training_accuracy = expand(rules.train_models_scale_after_split_10_iterations.output.training_accuracy, **config),
            test_accuracy = expand(rules.test_models_scale_after_split_10_iterations.output.test_accuracy, **config)
        output:
            scatter = os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_10_iterations.svg')
        conda:
            "../envs/python.yaml"
        params:
            colors = config['colors_complex']['model_colors']
        script:
            "../scripts/python/graphs/scatter_accuracies_10_iterations.py"

rule scatter_accuracies_10_iterations_only_hamming:
    """ Generate a connected scatter plot for each model to show their accuracy
        of prediction on the CV, training and test sets to highlight possible
        overfitting. Show the average accuracy plus std deviation across all 10
        iterations. These models were trained using only the combined_box_hamming
        feature."""
        input:
            cv_accuracy = expand(rules.hyperparameter_tuning_cv_scale_after_manual_split_only_hamming.output.best_hyperparameters, **config),
            training_accuracy = expand(rules.train_models_scale_after_manual_split_only_hamming.output.training_accuracy, **config),
            test_accuracy = expand(rules.test_models_scale_after_manual_split_only_hamming.output.test_accuracy, **config)
        output:
            scatter = os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_10_iterations_only_hamming.svg')
        conda:
            "../envs/python.yaml"
        params:
            colors = config['colors_complex']['model_colors']
        script:
            "../scripts/python/graphs/scatter_accuracies_10_iterations.py"

rule scatter_accuracies_20_iterations:
    """ Generate a connected scatter plot for each model to show their accuracy
        of prediction on the CV, training and test sets to highlight possible
        overfitting. Show the average accuracy plus std deviation across all 20
        iterations."""
        input:
            cv_accuracy = expand(os.path.join(config['path']['hyperparameter_tuning'],
                                    '{mod}_best_params_scale_after_split_{iteration_20}.tsv'),
                                    mod=['log_reg', 'svc', 'rf', 'gbm', 'knn'],
                                    iteration_20=config['iteration_20']),
            training_accuracy = expand(os.path.join(config['path']['training_accuracy'],
                                    '{mod}_training_accuracy_scale_after_split_{iteration_20}.tsv'),
                                    mod=['log_reg', 'svc', 'rf', 'gbm', 'knn'],
                                    iteration_20=config['iteration_20']),
            test_accuracy = expand(os.path.join(config['path']['test_accuracy'],
                                    '{mod}_test_accuracy_scale_after_split_{iteration_20}.tsv'),
                                    mod=['log_reg', 'svc', 'rf', 'gbm', 'knn'],
                                    iteration_20=config['iteration_20']),
        output:
            scatter = os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_20_iterations.svg')
        conda:
            "../envs/python.yaml"
        params:
            colors = config['colors_complex']['model_colors']
        script:
            "../scripts/python/graphs/scatter_accuracies_20_iterations.py"

rule scatter_accuracies_manual_split_iterations:
    """ Generate a connected scatter plot for each model to show their accuracy
        of prediction on the CV, training and test sets to highlight possible
        overfitting. Show the average accuracy plus std deviation across all 10
        manual split iterations."""
        input:
            cv_accuracy = expand(rules.hyperparameter_tuning_cv_scale_after_manual_split.output.best_hyperparameters, **config),
            training_accuracy = expand(rules.train_models_scale_after_manual_split.output.training_accuracy, **config),
            test_accuracy = expand(rules.test_models_scale_after_manual_split.output.test_accuracy, **config)
        output:
            scatter = os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_manual_split_iterations.svg')
        conda:
            "../envs/python.yaml"
        params:
            colors = config['colors_complex']['model_colors']
        script:
            "../scripts/python/graphs/scatter_accuracies_10_iterations.py"

rule roc_curve:
    """ Generate a roc curve for each trained model based on their performance
        on the test dataset."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        log_reg = os.path.join(config['path']['trained_models'],
                                    'log_reg_trained.sav'),
        svc = os.path.join(config['path']['trained_models'],
                                    'svc_trained.sav'),
        rf = os.path.join(config['path']['trained_models'],
                                    'rf_trained.sav'),
        gbm = os.path.join(config['path']['trained_models'],
                                    'gbm_trained.sav'),
        knn = os.path.join(config['path']['trained_models'],
                                    'knn_trained.sav'),
    output:
        roc_curve = os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models.svg')
    script:
        "../scripts/python/graphs/roc_curve.py"

rule roc_curve_wo_clusters:
    """ Generate a roc curve for each trained model based on their performance
        on the test dataset (CV, train and test sets without SNORD113, 114, 115
        and 116 snoRNAs)."""
    input:
        df = rules.remove_snoRNA_clusters.output.df_wo_clusters,
        log_reg = os.path.join(config['path']['trained_models'],
                                    'log_reg_trained_wo_clusters.sav'),
        svc = os.path.join(config['path']['trained_models'],
                                    'svc_trained_wo_clusters.sav'),
        rf = os.path.join(config['path']['trained_models'],
                                    'rf_trained_wo_clusters.sav'),
        gbm = os.path.join(config['path']['trained_models'],
                                    'gbm_trained_wo_clusters.sav'),
        knn = os.path.join(config['path']['trained_models'],
                                    'knn_trained_wo_clusters.sav'),
    output:
        roc_curve = os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_wo_clusters.svg')
    script:
        "../scripts/python/graphs/roc_curve_wo_clusters.py"

rule roc_curve_wo_feature:
    """ Generate a roc curve for each trained model based on their performance
        on the test dataset (CV, train and test without the top 3 most
        predictive features)."""
    input:
        df = rules.remove_feature.output.df_wo_feature,
        log_reg = os.path.join(config['path']['trained_models'],
                                    'log_reg_trained_wo_{feature_effect}.sav'),
        svc = os.path.join(config['path']['trained_models'],
                                    'svc_trained_wo_{feature_effect}.sav'),
        rf = os.path.join(config['path']['trained_models'],
                                    'rf_trained_wo_{feature_effect}.sav'),
        gbm = os.path.join(config['path']['trained_models'],
                                    'gbm_trained_wo_{feature_effect}.sav'),
        knn = os.path.join(config['path']['trained_models'],
                                    'knn_trained_wo_{feature_effect}.sav'),
    output:
        roc_curve = os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_wo_{feature_effect}.svg')
    script:
        "../scripts/python/graphs/roc_curve_wo_feature.py"

rule roc_curve_one_feature:
    """ Generate a roc curve for each trained model based on their performance
        on the test dataset (CV, train and test with only one feature (between
        one of the most predictive features)."""
    input:
        df = rules.keep_one_feature.output.df_one_feature,
        log_reg = os.path.join(config['path']['trained_models'],
                                    'log_reg_trained_only_{one_feature}.sav'),
        svc = os.path.join(config['path']['trained_models'],
                                    'svc_trained_only_{one_feature}.sav'),
        rf = os.path.join(config['path']['trained_models'],
                                    'rf_trained_only_{one_feature}.sav'),
        gbm = os.path.join(config['path']['trained_models'],
                                    'gbm_trained_only_{one_feature}.sav'),
        knn = os.path.join(config['path']['trained_models'],
                                    'knn_trained_only_{one_feature}.sav'),
    output:
        roc_curve = os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_only_{one_feature}.svg')
    script:
        "../scripts/python/graphs/roc_curve_one_feature.py"

rule roc_curve_scale_after_split:
    """ Generate a roc curve for each trained model based on their performance
        on the test dataset based on the CV, train, test sets scaled after they
        were generated."""
    input:
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        y_test = rules.fill_na_feature_scaling_after_split.output.y_test,
        log_reg = os.path.join(config['path']['trained_models'],
                                    'log_reg_trained_scale_after_split.sav'),
        svc = os.path.join(config['path']['trained_models'],
                                    'svc_trained_scale_after_split.sav'),
        rf = os.path.join(config['path']['trained_models'],
                                    'rf_trained_scale_after_split.sav'),
        gbm = os.path.join(config['path']['trained_models'],
                                    'gbm_trained_scale_after_split.sav'),
        knn = os.path.join(config['path']['trained_models'],
                                    'knn_trained_scale_after_split.sav'),
    output:
        roc_curve = os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_split.svg')
    script:
        "../scripts/python/graphs/roc_curve_scale_after_split.py"

rule roc_curve_scale_after_split_10_iterations:
    """ Generate a roc curve for each trained model based on their performance
        on the test dataset based on the CV, train, test sets scaled after they
        were generated. Add a shadow of +/- 1 stdev around each model roc curve."""
    input:
        X_test = expand(rules.fill_na_feature_scaling_after_split_10_iterations.output.test, **config),
        y_test = expand(rules.fill_na_feature_scaling_after_split_10_iterations.output.y_test, **config),
        pickled_trained_model = expand(rules.train_models_scale_after_split_10_iterations.output.pickled_trained_model,
                                    iteration=config['iteration'], models2=config['models2'])
    output:
        roc_curve = os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_split_10_iterations.svg')
    params:
        model_colors_dict = config['colors_complex']['model_colors']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/roc_curve_scale_after_split_10_iterations.py"

rule roc_curve_scale_after_manual_split:
    """ Generate a roc curve for each trained model based on their performance
        on the test dataset based on the CV, train, test sets scaled after they
        were generated. Add a shadow of +/- 1 stdev around each model roc curve."""
    input:
        X_test = expand(rules.fill_na_feature_scaling_after_manual_split.output.test, **config),
        y_test = expand(rules.fill_na_feature_scaling_after_manual_split.output.y_test, **config),
        pickled_trained_model = expand(rules.train_models_scale_after_manual_split.output.pickled_trained_model,
                                    manual_iteration=config['manual_iteration'], models2=config['models2'])
    output:
        roc_curve = os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_manual_split.svg')
    params:
        model_colors_dict = config['colors_complex']['model_colors']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/roc_curve_scale_after_split_10_iterations.py"

rule roc_curve_scale_after_split_10_iterations_only_hamming:
    """ Generate a roc curve for each trained model based on their performance
        on the test dataset based on the CV, train, test sets scaled after they
        were generated. Add a shadow of +/- 1 stdev around each model roc curve.
        These models were trained only with the feature combined_box_hamming."""
    input:
        X_test = expand(rules.fill_na_feature_scaling_after_manual_split_only_hamming.output.test, **config),
        y_test = expand(rules.fill_na_feature_scaling_after_manual_split_only_hamming.output.y_test, **config),
        pickled_trained_model = expand(rules.train_models_scale_after_manual_split_only_hamming.output.pickled_trained_model,
                                    manual_iteration=config['manual_iteration'], models2=config['models2'])
    output:
        roc_curve = os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_split_10_iterations_only_hamming.svg')
    params:
        model_colors_dict = config['colors_complex']['model_colors']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/roc_curve_scale_after_split_10_iterations.py"

rule shap_global:
    """ Generate a summary plot based on SHAP values (Shapley additive
        explanations) for each model explaining the effect of all features on
        the prediction of snoRNA abundance status in the test set."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    output:
        #summary_plot = os.path.join(config['figures']['summary_shap'],
        #                '{models}_all_features_test_set.svg'),
        summary_plot = os.path.join(config['figures']['summary_shap'],
                        '{models}_all_features_test_set_100_background.svg')
    script:
        "../scripts/python/graphs/summary_shap.py"

rule shap_global_snotype:
    """ Generate a summary plot based on SHAP values (Shapley additive
        explanations) for each model explaining the effect of all features on
        the prediction of either C/D or H/ACA snoRNA abundance status in the
        test set."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    output:
        summary_plot = os.path.join(config['figures']['summary_shap_snotype'],
                        '{models}_{sno_type}_all_features_test_set_100_background.svg')
    script:
        "../scripts/python/graphs/summary_shap_snotype.py"

rule shap_global_snotype_scale_after_split:
    """ Generate a summary plot based on SHAP values (Shapley additive
        explanations) for each model explaining the effect of all features on
        the prediction of either C/D or H/ACA snoRNA abundance status in the
        test set."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        df = config['path']['feature_df'],
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    output:
        summary_plot = os.path.join(config['figures']['summary_shap_snotype'],
                        '{models2}_{sno_type}_test_set_100_background_scale_after_split.svg')
    script:
        "../scripts/python/graphs/summary_shap_snotype_scale_after_split.py"

rule shap_global_hg_biotype:
    """ Generate a summary plot based on SHAP values (Shapley additive
        explanations) for each model explaining the effect of all features on
        the prediction of intronic or intergenic snoRNA abundance status in the
        test set."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    output:
        summary_plot = os.path.join(config['figures']['summary_shap_hg_biotype'],
                        '{models}_{hg_biotype}_all_features_test_set_100_background.svg')
    script:
        "../scripts/python/graphs/summary_shap_hg_biotype.py"

rule shap_global_hg_biotype_scale_after_split:
    """ Generate a summary plot based on SHAP values (Shapley additive
        explanations) for each model explaining the effect of all features on
        the prediction of intronic or intergenic snoRNA abundance status in the
        test set."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        df = config['path']['feature_df'],
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    output:
        summary_plot = os.path.join(config['figures']['summary_shap_hg_biotype'],
                        '{models2}_{hg_biotype}_test_set_100_background_scale_after_split.svg')
    script:
        "../scripts/python/graphs/summary_shap_hg_biotype_scale_after_split.py"

rule shap_local:
    """ Explain locally (per snoRNA) the decision taken by the classifier using
        SHAP on the examples in the test set."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    output:
        decision_plot = os.path.join(config['figures']['decision_plot'],
                        'ENSG00000212498_{models}_all_features_test_set_100_background.svg')
    script:
        "../scripts/python/graphs/decision_plot.py"

rule global_bar_plot:
    """ Create a bar chart of average feature importance (across all snoRNAs in
        test set) based on SHAP values per classifier"""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    output:
        bar_plot = os.path.join(config['figures']['bar'],
                        '{models}_global_feature_importance.svg')
    script:
        "../scripts/python/graphs/global_shap_bar_plot.py"

rule global_bar_plot_scale_after_split:
    """ Create a bar chart of average feature importance (across all snoRNAs in
        test set) based on SHAP values per classifier (dataset scaled after the split)"""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    output:
        bar_plot = os.path.join(config['figures']['bar'],
                        '{models2}_global_feature_importance_scale_after_split.svg')
    script:
        "../scripts/python/graphs/global_shap_bar_plot_scale_after_split.py"

rule upset_models_confusion:
    """ Create an upset plot per confusion matrix value (TN, FN, FP and TP) to
        compare the number of common snoRNAs (mis)classified by each model.
        Return also a merged df of all confusion matrices values (for all models
        in one df) (output hardcoded within script)."""
    input:
        log_reg = os.path.join(config['path']['confusion_matrix_f1'],
                        'log_reg_confusion_matrix_values_per_sno.tsv'),
        svc = os.path.join(config['path']['confusion_matrix_f1'],
                        'svc_confusion_matrix_values_per_sno.tsv'),
        rf = os.path.join(config['path']['confusion_matrix_f1'],
                        'rf_confusion_matrix_values_per_sno.tsv'),
        gbm = os.path.join(config['path']['confusion_matrix_f1'],
                        'gbm_confusion_matrix_values_per_sno.tsv')
    output:
        upset = os.path.join(config['figures']['upset'],
                        '{confusion_value}_all_models.svg'),
    script:
        "../scripts/python/graphs/upset_per_confusion_value.py"

rule upset_models_confusion_scale_after_split:
    """ Create an upset plot per confusion matrix value (TN, FN, FP and TP) to
        compare the number of common snoRNAs (mis)classified by each model.
        Return also a merged df of all confusion matrices values (for all models
        in one df) (output hardcoded within script). RF is excluded due to overfitting."""
    input:
        log_reg = os.path.join(config['path']['confusion_matrix_f1'],
                        'log_reg_confusion_matrix_values_per_sno_scale_after_split.tsv'),
        svc = os.path.join(config['path']['confusion_matrix_f1'],
                        'svc_confusion_matrix_values_per_sno_scale_after_split.tsv'),
        gbm = os.path.join(config['path']['confusion_matrix_f1'],
                        'gbm_confusion_matrix_values_per_sno_scale_after_split.tsv'),
        knn = os.path.join(config['path']['confusion_matrix_f1'],
                        'knn_confusion_matrix_values_per_sno_scale_after_split.tsv')
    output:
        upset = os.path.join(config['figures']['upset'],
                        '{confusion_value}_all_models_scale_after_split.svg')
    script:
        "../scripts/python/graphs/upset_per_confusion_value.py"

rule upset_models_confusion_scale_after_split_10_iterations:
    """ Create an upset plot per confusion matrix value (TN, FN, FP and TP) to
        compare the number of common snoRNAs (mis)classified by each model.
        Return also a merged df of all confusion matrices values (for all models
        in one df) (output hardcoded within script). GBM and KNN are excluded due to overfitting."""
    input:
        log_reg = expand(rules.confusion_matrix_f1_scale_after_split_10_iterations.output.info_df,
                        iteration=config['iteration'], models2='log_reg', allow_missing=True),
        svc = expand(rules.confusion_matrix_f1_scale_after_split_10_iterations.output.info_df,
                        iteration=config['iteration'], models2='svc', allow_missing=True),
        rf = expand(rules.confusion_matrix_f1_scale_after_split_10_iterations.output.info_df,
                        iteration=config['iteration'], models2='rf', allow_missing=True)
    params:
        df_output_path = expand("results/tables/confusion_matrix_f1/merged_confusion_matrix_{iteration}.tsv",
                                iteration=config['iteration'])
    output:
        upset = os.path.join(config['figures']['upset'],
                        '{confusion_value}_all_models_scale_after_split_10_iterations.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/upset_per_confusion_value_10_iterations.py"

rule upset_models_confusion_scale_after_manual_split:
    """ Create an upset plot per confusion matrix value (TN, FN, FP and TP) to
        compare the number of common snoRNAs (mis)classified by each model.
        Return also a merged df of all confusion matrices values (for all models
        in one df) (output hardcoded within script). GBM and KNN are excluded due to overfitting."""
    input:
        log_reg = expand(rules.confusion_matrix_f1_scale_after_manual_split.output.info_df,
                        manual_iteration=config['manual_iteration'], models2='log_reg', allow_missing=True),
        svc = expand(rules.confusion_matrix_f1_scale_after_manual_split.output.info_df,
                        manual_iteration=config['manual_iteration'], models2='svc', allow_missing=True),
        rf = expand(rules.confusion_matrix_f1_scale_after_manual_split.output.info_df,
                        manual_iteration=config['manual_iteration'], models2='rf', allow_missing=True)
    params:
        df_output_path = expand("results/tables/confusion_matrix_f1/merged_confusion_matrix_{manual_iteration}.tsv",
                                manual_iteration=config['manual_iteration'])
    output:
        upset = os.path.join(config['figures']['upset'],
                        '{confusion_value}_all_models_scale_after_manual_split.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/upset_per_confusion_value_10_iterations.py"

rule upset_top_features_df:
    """ Create a dataframe containing the top 5 features (their name, rank and
        feature_importance) of all four models (not RF because it overfits)."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        log_reg = os.path.join(config['path']['trained_models'],
                                    'log_reg_trained_scale_after_split.sav'),
        other_model = expand(os.path.join(config['path']['trained_models'],
                                    '{mod}_trained_scale_after_split.sav'),
                                    mod = ['svc', 'gbm', 'knn']),
    output:
        top_features_df = config['path']['top_5_features_df']
    script:
        "../scripts/python/graphs/upset_top_features_df.py"

rule all_feature_rank_df:
    """ Create a dataframe containing the rank of importance for each feature
        and per model."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        log_reg = os.path.join(config['path']['trained_models'],
                                    'log_reg_trained_scale_after_split.sav'),
        other_model = expand(os.path.join(config['path']['trained_models'],
                                    '{mod}_trained_scale_after_split.sav'),
                                    mod = ['svc', 'gbm', 'knn'])
    output:
        rank_features_df = config['path']['all_feature_rank_df']
    script:
        "../scripts/python/all_feature_rank_df.py"

rule all_feature_rank_df_10_iterations:
    """ Create a dataframe containing the rank of importance for each feature
        and per model for all 10 iterations."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split_10_iterations.output.train,
        X_test = rules.fill_na_feature_scaling_after_split_10_iterations.output.test,
        log_reg = expand(rules.train_models_scale_after_split_10_iterations.output.pickled_trained_model,
                                    models2='log_reg', allow_missing=True),
        svc = expand(rules.train_models_scale_after_split_10_iterations.output.pickled_trained_model,
                                    models2='svc', allow_missing=True),
        rf = expand(rules.train_models_scale_after_split_10_iterations.output.pickled_trained_model,
                                    models2='rf', allow_missing=True)
    output:
        rank_features_df = "results/tables/all_features_rank_across_models_{iteration}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/all_feature_rank_df_10_iterations.py"

rule concat_feature_rank_iterations_df:
    """ Concat all dfs produced by the all_feature_rank_df_10_iterations rule
        into one df containing all iterations."""
    input:
        dfs = expand(rules.all_feature_rank_df_10_iterations.output.rank_features_df,
                    iteration=config['iteration'])
    output:
        concat_df = config['path']['all_feature_rank_df_10_iterations']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/concat_feature_rank_iterations_df.py"


rule violin_feature_rank:
    """ Create a violin plot where each violin on the x axis corresponds to one
        feature and where the distribution of ranks (across 4 models) is
        represented on the y axis."""
    input:
        rank_features_df = rules.all_feature_rank_df.output.rank_features_df
    output:
        violin = os.path.join(config['figures']['violin'], 'ranks_per_feature.svg')
    params:
        model_colors = config['colors_complex']['model_colors']
    script:
        "../scripts/python/graphs/violin_feature_rank.py"

rule violin_feature_rank_10_iterations:
    """ Create a violin plot where each violin on the x axis corresponds to one
        feature and where the distribution of ranks (across 3 models and their
        10 respective iterations) is represented on the y axis."""
    input:
        rank_features_df = rules.concat_feature_rank_iterations_df.output.concat_df
    output:
        violin = os.path.join(config['figures']['violin'], 'ranks_per_feature_10_iterations.svg')
    params:
        model_colors = config['colors_complex']['model_colors']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/violin_feature_rank_10_iterations.py"

rule heatmap_feature_rank_correlation:
    """ Create a heatmap showing the spearman correlation between feature
        predictive ranks across 10 iterations. One heatmap per model."""
    input:
        rank_features_df = rules.concat_feature_rank_iterations_df.output.concat_df
    output:
        heatmaps = expand(os.path.join(config['figures']['heatmap'],
                    '{model_name}_feature_rank_correlation_10_iterations.svg'),
                    model_name=['log_reg', 'svc', 'rf'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/heatmap_feature_rank_correlation.py"

rule clustermap_feature_rank_correlation:
    """ Create a clustered heatmap showing the spearman correlation between feature
        predictive ranks across 10 iterations. One heatmap per model."""
    input:
        rank_features_df = rules.concat_feature_rank_iterations_df.output.concat_df
    output:
        heatmaps = expand(os.path.join(config['figures']['heatmap'],
                    '{model_name}_feature_rank_correlation_10_iterations_clustered.svg'),
                    model_name=['log_reg', 'svc', 'rf'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/clustermap_feature_rank_correlation.py"

rule pairplot_top_5:
    """ Create a pairplot per snoRNA type to show the dependence between the top
        5 features (across models)."""
    input:
        rank_features_df = rules.all_feature_rank_df.output.rank_features_df,
        all_features_df = config['path']['one_hot_encoded_df']
    output:
        pairplot_cd = os.path.join(config['figures']['pairplot'], 'top_5_features_cd.svg'),
        pairplot_haca = os.path.join(config['figures']['pairplot'], 'top_5_features_haca.svg')
    params:
        hue_color = config['colors_complex']['label']
    script:
        "../scripts/python/graphs/pairplot_top_5.py"

rule upset_top_features:
    """ Create with UpSetR an upset plot of the intersection between the top 5
        most predictive features of all models (without RF because it overfits).
        It takes the df of upset_top_features_df as input."""
    input:
        df = rules.upset_top_features_df.output.top_features_df
    output:
        upset = os.path.join(config['figures']['upset'],
                        'top_5_features_intersection_all_models_scale_after_split.svg')
    conda:
        "../envs/upset_r.yaml"
    script:
        "../scripts/r/upset_top_features.R"

rule donut_top_features_percent:
    """ Create a donut chart for each intersection category (feature common to 4
        models, 3, 2 or 1 model) of top 5 features across all four models. The
        percentage represents the proportion of feature(s) per intersection
        category that are ranked the 1st, 2nd, ..., 5th most predictive feature
        across models."""
    input:
        df = rules.upset_top_features_df.output.top_features_df
    output:
        donut = os.path.join(config['figures']['donut'],
                        'top_feature_intersection_rank_percent.svg')
    params:
        rank_colors = config['colors_complex']['rank_top_features']
    script:
        "../scripts/python/graphs/donut_top_features_rank_percent.py"

rule shap_local_FN:
    """ Explain locally (per snoRNA) the decision taken by the 4 classifiers
        using SHAP on the false negatives (FN) examples common to all models.
        The graphs are created within the script so that we don't have to create
        a time-consuming explainer for all FN snoRNAs (figure output not
        explicited here in the output)."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    params:
        false_negatives = config['FN'],
        decision_plot_FN = config['figures']['decision_plot_FN']
    output:
        shap_local_FN_log = "log/shap_local_FN_{models}.log"
    script:
        "../scripts/python/graphs/decision_plot_FN.py"

rule shap_local_FP:
    """ Explain locally (per snoRNA) the decision taken by the 4 classifiers
        using SHAP on the false positives (FP) examples common to all models.
        The graphs are created within the script so that we don't have to create
        a time-consuming explainer for all FP snoRNAs (figure output not
        explicited here in the output)."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    params:
        false_positives = config['FP'],
        decision_plot_FP = config['figures']['decision_plot_FP']
    output:
        shap_local_FP_log = "log/shap_local_FP_{models}.log"
    script:
        "../scripts/python/graphs/decision_plot_FP.py"

rule shap_local_TN:
    """ Explain locally (per snoRNA) the decision taken by the 4 classifiers
        using SHAP on the true negatives (TN) examples common to all models.
        The graphs are created within the script so that we don't have to create
        a time-consuming explainer for all TN snoRNAs (figure output not
        explicited here in the output)."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    params:
        true_negatives = config['TN'],
        decision_plot_TN = config['figures']['decision_plot_TN']
    output:
        shap_local_TN_log = "log/shap_local_TN_{models}.log"
    script:
        "../scripts/python/graphs/decision_plot_TN.py"

rule shap_local_TP:
    """ Explain locally (per snoRNA) the decision taken by the 4 classifiers
        using SHAP on the true positives (TP) examples common to all models.
        The graphs are created within the script so that we don't have to create
        a time-consuming explainer for all TP snoRNAs (figure output not
        explicited here in the output)."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    params:
        true_positives = config['TP'],
        decision_plot_TP = config['figures']['decision_plot_TP']
    output:
        shap_local_TP_log = "log/shap_local_TP_{models}.log"
    script:
        "../scripts/python/graphs/decision_plot_TP.py"

rule shap_local_FN_scale_after_split:
    """ Explain locally (per snoRNA) the decision taken by the 4 classifiers
        using SHAP on the false negatives (FN) examples common to all models.
        The graphs are created within the script so that we don't have to create
        a time-consuming explainer for all FN snoRNAs (figure output not
        explicited here in the output)."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    params:
        false_negatives = config['FN'],
        decision_plot_FN = config['figures']['decision_plot_FN_scale_after_split']
    output:
        shap_local_FN_log = "log/shap_local_FN_{models2}_scale_after_split.log"
    script:
        "../scripts/python/graphs/decision_plot_FN_scale_after_split.py"

rule shap_local_FP_scale_after_split:
    """ Explain locally (per snoRNA) the decision taken by the 4 classifiers
        using SHAP on the false positives (FP) examples common to all models.
        The graphs are created within the script so that we don't have to create
        a time-consuming explainer for all FP snoRNAs (figure output not
        explicited here in the output)."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    params:
        false_positives = config['FP'],
        decision_plot_FP = config['figures']['decision_plot_FP_scale_after_split']
    output:
        shap_local_FP_log = "log/shap_local_FP_{models2}_scale_after_split.log"
    script:
        "../scripts/python/graphs/decision_plot_FP_scale_after_split.py"

rule shap_local_TN_scale_after_split:
    """ Explain locally (per snoRNA) the decision taken by the 4 classifiers
        using SHAP on the true negatives (TN) examples common to all models.
        The graphs are created within the script so that we don't have to create
        a time-consuming explainer for all TN snoRNAs (figure output not
        explicited here in the output)."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    params:
        true_negatives = config['TN'],
        decision_plot_TN = config['figures']['decision_plot_TN_scale_after_split']
    output:
        shap_local_TN_log = "log/shap_local_TN_{models2}_scale_after_split.log"
    script:
        "../scripts/python/graphs/decision_plot_TN_scale_after_split.py"

rule shap_local_TP_scale_after_split:
    """ Explain locally (per snoRNA) the decision taken by the 4 classifiers
        using SHAP on the true positives (TP) examples common to all models.
        The graphs are created within the script so that we don't have to create
        a time-consuming explainer for all TP snoRNAs (figure output not
        explicited here in the output)."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    params:
        true_positives = config['TP'],
        decision_plot_TP = config['figures']['decision_plot_TP_scale_after_split']
    output:
        shap_local_TP_log = "log/shap_local_TP_{models2}_scale_after_split.log"
    script:
        "../scripts/python/graphs/decision_plot_TP_scale_after_split.py"


rule structure_snora77b:
    """ Generate with RNAfold the structure of SNORA77B (ENSG00000264346) and
        its supposed terminal stem using the snoRNA sequence and its flanking
        intronic sequences combined as one fasta file."""
    input:
        sno_sequences = rules.format_snoRNA_sequences.output.sno_sequences,
        flanking = rules.get_fasta.output.sequences
    params:
        bash_script = "./scripts/bash/snora77b_structure.sh"
    output:
        snora77b_terminal_stem = config['path']['snora77b_terminal_stem_fasta']
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "{params.bash_script} {input.sno_sequences} "
        "{input.flanking} {output.snora77b_terminal_stem} "

rule forgi_structure_snora77b_and_stem:
    """ Generate with forgi the structure of SNORA77B (ENSG00000264346) alone
        and also with its supposed terminal stem using the snoRNA sequence and
        its flanking intronic sequences."""
    input:
        snora77b_terminal_stem = config['path']['snora77b_terminal_stem_fasta'],
        all_snorna_structure = config['path']['structure_stability_fasta']
    output:
        snora77b_terminal_stem_figure = os.path.join(config['figures']['forgi_structure'],
                                "SNORA77B_terminal_stem.svg"),
        snora77b_figure = os.path.join(config['figures']['forgi_structure'],
                                "SNORA77B.svg"),
        snora77b_dot_bracket = config['path']['snora77b_fasta']
    script:
        "../scripts/python/graphs/forgi_snora77b.py"

rule bivariate_density_haca:
    """ Generate a bivariate density plot showing the difference between
        expressed and not expressed H/ACA snoRNAs according to their
        conservation and structure stability (MFE)."""
    input:
        all_features = config['path']['feature_df']
    params:
        ab_status_color = config['colors_complex']['abundance_cutoff_2']
    output:
        bivariate_density = os.path.join(config['figures']['bivariate_density'],
                                "HACA_abundance_status_conservation_mfe.svg")
    script:
        "../scripts/python/graphs/bivariate_density_haca.py"

rule violin_tpm_sno_type:
    """ Generate a violin plot of the log10(avg TPM across tissues) of snoRNAs
        separated by snoRNA type to highlight SNORA77B high abundance."""
    input:
        tpm_df = config['path']['sno_tpm_df_cutoff'],
        all_features = config['path']['feature_df']
    params:
        colors = config['colors_complex']['sno_type']
    output:
        violin = os.path.join(config['figures']['violin'], 'tpm_comparison_sno_type.svg')
    script:
        "../scripts/python/graphs/comparison_tpm_sno_type.py"


rule stem_comparison_tpm:
    """ Compare the abundance (TPM) of SNORA77B with other expressed H/ACA that
        have (or not) a good terminal stem stability. Compare also between C/D
        box snoRNAs. Create a violin plot out of it"""
    input:
        tpm_df = config['path']['sno_tpm_df_cutoff'],
        all_features = config['path']['feature_df']
    params:
        cd_colors = config['colors_complex']['terminal_stem_strength_cd'],
        haca_colors = config['colors_complex']['terminal_stem_strength_haca']
    output:
        violin_cd = os.path.join(config['figures']['violin'], 'cd_stem_comparison_tpm.svg'),
        violin_haca = os.path.join(config['figures']['violin'], 'haca_stem_comparison_tpm.svg')
    script:
        "../scripts/python/graphs/stem_comparison_tpm.py"

rule pca:
    """ Generate a PCA plot for the test dataset and for the whole initial
        dataset to see how much of the variance is explained by the first 2
        principal components."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df
    output:
        pca_all = os.path.join(config['figures']['pca'], "all_snoRNAs_{pca_hue}.svg"),
        pca_test = os.path.join(config['figures']['pca'], "snoRNAs_test_{pca_hue}.svg")
    params:
        colors_dict = lambda wildcards: config['colors_complex'][wildcards.pca_hue]
    script:
        "../scripts/python/graphs/PCA.py"

rule t_sne:
    """ Generate a t-SNE plot for the test dataset and for the whole initial
        dataset to see how much of the variance is explained by the components
        found by the tSNE."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df
    output:
        t_sne_all = os.path.join(config['figures']['t_sne'], "all_snoRNAs_{pca_hue}.svg"),
        t_sne_test = os.path.join(config['figures']['t_sne'], "snoRNAs_test_{pca_hue}.svg")
    params:
        colors_dict = lambda wildcards: config['colors_complex'][wildcards.pca_hue]
    script:
        "../scripts/python/graphs/tSNE.py"

rule global_feature_importance:
    """ Show global feature importance (implemented in sklearn for some models,
        not all) with a bar plot. Compare to SHAP global values."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained.sav')
    output:
        bar_plot = os.path.join(config['figures']['bar'],
                                    '{models}_sklearn_feature_importance.svg')
    script:
        "../scripts/python/graphs/feature_importance.py"

rule global_feature_importance_wo_top_10_all:
    """ Show global feature importance (implemented in sklearn for some models,
        not all) with a bar plot (without top_10_all features)"""
    input:
        df = os.path.join(config['path']['one_hot_encoded_df_wo_feature'],
                            'top_10_all.tsv'),
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained_wo_top_10_all.sav')
    output:
        bar_plot = os.path.join(config['figures']['bar'],
                                    '{models}_sklearn_feature_importance_wo_top_10_all.svg')
    script:
        "../scripts/python/graphs/feature_importance.py"

rule global_feature_importance_wo_top_11_all:
    """ Show global feature importance (implemented in sklearn for some models,
        not all) with a bar plot (without top_11_all features)"""
    input:
        df = os.path.join(config['path']['one_hot_encoded_df_wo_feature'],
                            'top_11_all.tsv'),
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained_wo_top_11_all.sav')
    output:
        bar_plot = os.path.join(config['figures']['bar'],
                                    '{models}_sklearn_feature_importance_wo_top_11_all.svg')
    script:
        "../scripts/python/graphs/feature_importance.py"

rule heatmap_shap:
    """ Create a heatmap of all snoRNAs in the test set (on the x axis)
        clustered by their explanation similarity, resulting in samples that
        have the same model output for the same reasons getting grouped together.
        On the y axis, you have the input features used by the models ordered by
        feature importance (also shown on the rigth side is the bar charts). The
        output of the model (f(x)) is shown above the heatmap."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        y_test = rules.fill_na_feature_scaling_after_split.output.y_test,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    params:
        labels_dict = config['colors_complex']['label']
    output:
        heatmap = os.path.join(config['figures']['heatmap'], '{models2}_shap_heatmap.svg')
    script:
        "../scripts/python/graphs/heatmap_shap.py"

rule sno_presence_in_all_test_sets:
    """ Find the snoRNAs that found in at leat 1 test set (across the 10
        iterations) and those that are never found in the test set."""
    input:
        test_sets = expand(rules.fill_na_feature_scaling_after_split_10_iterations.output.test, **config),
        all_sno_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df
    output:
        sno_presence_test_sets = config['path']['sno_presence_test_sets']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/sno_presence_in_all_test_sets.py"

rule regroup_sno_confusion_value_iterations:
    """ Regroup all snoRNAs by confusion value (TP, TN, FP, FN) across the 10
        iterations."""
    input:
        confusion_value_df = expand(rules.confusion_matrix_f1_scale_after_split_10_iterations.output.info_df,
                                models2=config['models3'], iteration=config['iteration'])
    output:
        sno_per_confusion_value = os.path.join(config['path']['sno_per_confusion_value'],
                                    '{confusion_value}_snoRNAs_10_iterations.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/regroup_sno_confusion_value_iterations.py"

rule num_feature_distribution_comparison_confusion_value_top10:
    """ Create feature distribution (density plot) to compare the numerical
        features contained in the top 10 most predictive features. """
    input:
        feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        sno_per_confusion_value = expand(rules.regroup_sno_confusion_value_iterations.output.sno_per_confusion_value,
                                            **config)
    output:
        density = os.path.join(config['figures']['density_confusion_value'],
                                "{comparison_confusion_val}_{top_10_numerical_features}.svg")
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/num_feature_distribution_comparison_confusion_value_top10.py"

rule num_feature_distribution_comparison_confusion_value_top10_snotype:
    """ Create feature distribution (density plot) to compare the numerical
        features contained in the top 10 most predictive features per snoRNA type. """
    input:
        feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        sno_per_confusion_value = expand(rules.regroup_sno_confusion_value_iterations.output.sno_per_confusion_value,
                                            **config)
    output:
        density = os.path.join(config['figures']['density_confusion_value'],
                                "{comparison_confusion_val}_{top_10_numerical_features}_{sno_type}.svg")
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/num_feature_distribution_comparison_confusion_value_top10_snotype.py"

rule cat_feature_distribution_comparison_confusion_value_top10:
    """ Create feature distribution (bar chart) to compare the categorical
        features contained in the top 10 most predictive features. """
    input:
        feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        sno_per_confusion_value = expand(rules.regroup_sno_confusion_value_iterations.output.sno_per_confusion_value,
                                            **config)
    output:
        bar = os.path.join(config['figures']['bar_confusion_value'],
                                "{comparison_confusion_val}_{top_10_categorical_features}.svg")
    params:
        color_dict = config['colors_complex']['top_10_categorical_features']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/cat_feature_distribution_comparison_confusion_value_top10.py"

rule cat_feature_distribution_comparison_confusion_value_top10_snotype:
    """ Create feature distribution (bar chart) to compare the categorical
        features contained in the top 10 most predictive features per snoRNA type. """
    input:
        feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        sno_per_confusion_value = expand(rules.regroup_sno_confusion_value_iterations.output.sno_per_confusion_value,
                                            **config)
    output:
        bar = os.path.join(config['figures']['bar_confusion_value'],
                                "{comparison_confusion_val}_{top_10_categorical_features}_{sno_type}.svg")
    params:
        color_dict = config['colors_complex']['top_10_categorical_features']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/cat_feature_distribution_comparison_confusion_value_top10_snotype.py"

rule get_all_shap_values:
    """ Get the SHAP values of all features for each snoRNA across the 10 test
        set iterations and for each model (log_reg, svc and rf). Shap values are
        given in the form of log(odds) not probability."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split_10_iterations.output.train,
        X_test = rules.fill_na_feature_scaling_after_split_10_iterations.output.test,
        pickled_trained_model = rules.train_models_scale_after_split_10_iterations.output.pickled_trained_model
    output:
        shap = os.path.join(config['path']['shap_10_iterations'], '{models2}_{iteration}_shap_values.tsv'),
        expected_value = os.path.join(config['path']['shap_10_iterations'], '{models2}_{iteration}_expected_value.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_all_shap_values.py"

rule get_all_shap_values_manual_split:
    """ Get the SHAP values of all features for each snoRNA across the 10 test
        set manual split iterations and for each model (log_reg, svc and rf). Shap values are
        given in the form of log(odds) not probability."""
    input:
        X_train = rules.fill_na_feature_scaling_after_manual_split.output.train,
        X_test = rules.fill_na_feature_scaling_after_manual_split.output.test,
        pickled_trained_model = rules.train_models_scale_after_manual_split.output.pickled_trained_model
    output:
        shap = os.path.join(config['path']['shap_10_iterations'], '{models2}_{manual_iteration}_shap_values.tsv'),
        expected_value = os.path.join(config['path']['shap_10_iterations'], '{models2}_{manual_iteration}_expected_value.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_all_shap_values_manual_split.py"

rule all_feature_rank_df_manual_split:
    """ Create a dataframe containing the rank of importance for each feature
        and per model for all 10 manual iterations."""
    input:
        shap_vals = expand(rules.get_all_shap_values_manual_split.output.shap,
                                    models2=config['models3'], allow_missing=True)
    output:
        rank_features_df = "results/tables/all_features_rank_across_models_{manual_iteration}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/all_feature_rank_df_manual_split.py"


rule concat_feature_rank_manual_split_iterations_df:
    """ Concat all dfs produced by the all_feature_rank_df_manual_split rule
        into one df containing all iterations."""
    input:
        dfs = expand(rules.all_feature_rank_df_manual_split.output.rank_features_df,
                    manual_iteration=config['manual_iteration'])
    output:
        concat_df = config['path']['all_feature_rank_df_manual_split']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/concat_feature_rank_iterations_df.py"

rule violin_feature_rank_manual_split:
    """ Create a violin plot where each violin on the x axis corresponds to one
        feature and where the distribution of ranks (across 3 models and their
        10 respective manual iterations) is represented on the y axis."""
    input:
        rank_features_df = rules.concat_feature_rank_manual_split_iterations_df.output.concat_df
    output:
        violin = os.path.join(config['figures']['violin'], 'ranks_per_feature_manual_split_iterations.svg')
    params:
        model_colors = config['colors_complex']['model_colors']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/violin_feature_rank_manual_split.py"

rule shap_global_snotype_scale_after_manual_split:
    """ Generate a summary plot based on SHAP values (Shapley additive
        explanations) for each model explaining the effect of all features on
        the prediction of either C/D or H/ACA snoRNA abundance status in the
        test set across 10 iterations combined."""
    input:
        fake_log = rules.modify_shap.output.fake_log,
        X_test = expand(rules.fill_na_feature_scaling_after_manual_split.output.test, **config),
        shap_values = expand(rules.get_all_shap_values_manual_split.output.shap,
                        manual_iteration=config['manual_iteration'], allow_missing=True),
        df = rules.merge_feature_df.output.feature_df
    output:
        summary_plot = os.path.join(config['figures']['summary_shap_snotype'],
                        '{models2}_{sno_type}_test_set_scale_after_manual_split.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/summary_shap_snotype_scale_after_manual_split.py"

rule bar_shap_global_snotype_scale_after_split_10_iterations:
    """ Generate a bar plot based on SHAP values (Shapley additive
        explanations) for each model explaining the effect of all features on
        the prediction of either C/D or H/ACA snoRNA abundance status in the
        test set across 10 iterations combined."""
    input:
        fake_log = rules.modify_shap.output.fake_log,
        X_test = expand(rules.fill_na_feature_scaling_after_manual_split.output.test, **config),
        shap_values = expand(rules.get_all_shap_values_manual_split.output.shap,
                        manual_iteration=config['manual_iteration'], allow_missing=True),
        df = rules.merge_feature_df.output.feature_df
    output:
        summary_plot = os.path.join(config['figures']['bar'],
                        '{models2}_{sno_type}_test_set_scale_after_manual_split.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/bar_summary_shap_snotype_scale_after_manual_split.py"

rule decision_plot_scale_after_split_10_iterations_per_confusion_value:
    """ Explain locally (per snoRNA) the decision taken by the 3 classifiers
        using SHAP on the confusion value (FN, FP, TN or TP) examples present across the 10
        iterations and that are common to all models."""
    input:
        shap_values = expand(rules.get_all_shap_values.output.shap, models2=config['models3'], iteration=config['iteration']),
        expected_value = expand(rules.get_all_shap_values.output.expected_value, models2=config['models3'], iteration=config['iteration']),
        sno_per_confusion_value = expand(rules.regroup_sno_confusion_value_iterations.output.sno_per_confusion_value, **config)
    output:
        shap_local_log_reg = os.path.join(config['figures']['decision_plot_10_iterations'], '{confusion_value}_10_iterations_log_reg.svg'),
        shap_local_svc = os.path.join(config['figures']['decision_plot_10_iterations'], '{confusion_value}_10_iterations_svc.svg'),
        shap_local_rf = os.path.join(config['figures']['decision_plot_10_iterations'], '{confusion_value}_10_iterations_rf.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/decision_plot_scale_after_split_10_iterations_per_confusion_value.py"

rule concat_all_shap_values:
    """ Concat all the SHAP values across the 10 test set iterations into one
        per model and confusion value."""
    input:
        all_shap = expand(rules.get_all_shap_values.output.shap, iteration=config['iteration'], allow_missing=True),
        sno_per_confusion_value = expand(rules.regroup_sno_confusion_value_iterations.output.sno_per_confusion_value, **config)
    output:
        concat_shap = os.path.join(config['path']['shap_10_iterations'], '{models2}_{confusion_value}_shap_values_concat.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/concat_all_shap_values.py"

rule density_FP_vs_TN_host_expressed:
    """ Compare the distribution of terminal_stem_mfe, combined_box_hamming and
        sno_mfe between all real FP (which all have an expressed HG) and all TN
        that have an expressed HG. """
    input:
        sno_per_confusion_value = expand(rules.regroup_sno_confusion_value_iterations.output.sno_per_confusion_value, **config),
        all_features_df = config['path']['feature_df']
    output:
        density = os.path.join(config['figures']['density_confusion_value'],
                                "host_expressed_FP_vs_TN_{top_10_numerical_features}.svg")
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/density_FP_vs_TN_host_expressed.py"

rule real_confusion_value_df:
    """ For each confusion value, get only the real snoRNAs (those always
        predicted as such) and their feature values."""
    input:
        sno_per_confusion_value = expand(rules.regroup_sno_confusion_value_iterations.output.sno_per_confusion_value, **config),
        all_features_df = config['path']['feature_df']
    output:
        real_confusion_value_df = os.path.join(config['path']['real_confusion_value'], '{confusion_value}_w_features.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/real_confusion_value_df.py"

rule violin_tpm_FN_TP:
    """ Create a violin plot of the log2(TPM) for all real false negatives (FN)
        vs all real true positives (TP) vs all TN vs all FP ("real" meaning that these snoRNAs are
        always predicted as their confusion value across models and iterations)."""
    input:
        sno_per_confusion_value = expand(rules.real_confusion_value_df.output.real_confusion_value_df,
                                            **config),
        tpm_df = config['path']['sno_tpm_df_cutoff']
    output:
        violin = os.path.join(config['figures']['violin'], "FN_vs_TP_vs_TP_vs_FP_tpm.svg")
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/violin_tpm_FN_TP.py"
