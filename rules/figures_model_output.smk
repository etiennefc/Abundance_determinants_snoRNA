import os

include: "cv_train_test.smk"
include: "cv_train_test_wo_clusters.smk"
include: "cv_train_test_wo_feature.smk"
include: "cv_train_test_one_feature.smk"
include: "cv_train_test_scale_after_split.smk"
include: "terminal_stem.smk"
include: "structure.smk"

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
                                    mod = ['svc', 'gbm', 'knn']),
    output:
        rank_features_df = config['path']['all_feature_rank_df']
    script:
        "../scripts/python/all_feature_rank_df.py"

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
