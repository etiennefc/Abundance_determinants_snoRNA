#!/usr/bin/python3
from snakemake.io import expand
import os
def get_figures_path(config):
    """Return list of figures (and their path) to generate from config"""
    files = []
    numerical = "{numerical_features}_{feature_hue}.svg"
    cd_numerical = "numerical_features_cd.svg"
    haca_numerical = "numerical_features_haca.svg"

    # Expressed vs not expressed RNA per biotype
    #files.append(os.path.join(config['figures']['bar'],
    #            'abundance_status_per_biotype.svg'))

    # Connected scatter plot of average accuracies of all models in CV, train and test sets across 10 iterations
    #files.append(os.path.join(config['figures']['scatter'], 'all_model_accuracies_cv_train_test_10_iterations.svg'))
    #files.extend(expand(config['figures']['upset'] + '{confusion_value}_all_models_scale_after_split_10_iterations.svg', **config))

    # ROC curve and connected scatter plot of average accuracies of all models in CV, train and test sets across 10 iterations
    # trained only with combined_box_hamming
    files.append(os.path.join(config['figures']['roc'],
                'roc_curves_test_set_5_models_scale_after_split_10_iterations_only_hamming.svg'))
    files.append(os.path.join(config['figures']['scatter'], 'all_model_accuracies_cv_train_test_10_iterations_only_hamming.svg'))


    # ROC curve and connected scatter plot of average accuracies of all models in CV, train and test sets across 10 manual split iterations
    files.append(os.path.join(config['figures']['roc'],
                    'roc_curves_test_set_5_models_scale_after_manual_split.svg'))
    files.append(os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_manual_split_iterations.svg'))

    # Violin plot of feature rank across 10 manual split iterations and upset plot of confusion value
    files.append(os.path.join(config['figures']['violin'], 'ranks_per_feature_manual_split_iterations.svg'))
    files.extend(expand(os.path.join(config['figures']['upset'],
                    '{confusion_value}_all_models_scale_after_manual_split.svg'), **config))

    # Pie chart of the number of snoRNAs per confusion value (TP, TN, FP, FN)
    #files.append(os.path.join(config['figures']['pie'], 'sno_number_per_confusion_value.svg'))

    # Donut chart of the number of snoRNAs per confusion value (TP, TN, FP, FN) (outer donut) and their host biotype (inner donut)
    files.append(os.path.join(config['figures']['donut'],
                        'confusion_value_host_biotype.svg'))

    # SHAP summary plot of all 10 manual split iterations per sno type (dot plot and bar plot)
    files.extend(expand(os.path.join(config['figures']['summary_shap_snotype'],
                                '{models2}_{sno_type}_test_set_scale_after_manual_split.svg'), **config))
    files.extend(expand(os.path.join(config['figures']['bar'],
                    '{models2}_{sno_type}_test_set_scale_after_manual_split.svg'), **config))

    # Multi Decision plot of snoRNAs in GAS5 or SNHG17
    files.extend(expand(os.path.join(config['figures']['decision_plot'], '{multi_HG_diff_label}_snoRNAs_{models2}.svg'), **config))

    # Pie and donut charts for mono vs multi HG with different snoRNA labels within the same HG
    files.append(os.path.join(config['figures']['pie'], 'mono_vs_multi_HG.svg'))

    # Decision plots for potential interesting snoRNAS (ex: SNORA77B, SNORD86)
    files.extend(expand(os.path.join(config['figures']['decision_plot_interesting_snoRNAs'], '{interesting_sno_ids}_{models2}.svg'), **config))

    # Density plots to compare per confusion_value the top 10 numerical features distributions
    #files.extend(expand(os.path.join(config['figures']['density_confusion_value'],
    #                            "confusion_value_comparison_{top_10_numerical_features}.svg"), **config))

    # Bar plots to compare per confusion_value the top 10 categorical features distributions
    files.extend(expand(os.path.join(config['figures']['bar_confusion_value'],
                                "confusion_value_comparison_{top_10_categorical_features}.svg"), **config))

    # Violin plot to show the log2(avg TPM) per confusion value
    files.append(os.path.join(config['figures']['violin'], "avg_tpm_per_confusion_value.svg"))

    # Horizontal stacked bar chart to show the number of sno per HG (per confusion value)
    files.append(os.path.join(config['figures']['hbar'], 'nb_sno_per_confusion_val_multi_HG.svg'))

    # ROC curve of models trained with top3 or top4 features
    files.append(os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_manual_split_top3.svg'))
    files.append(os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_manual_split_top4.svg'))

    # Connected scatter plot of average accuracies of all models in CV, train and test sets across 10 iterations of models trained on top3 of top4 features
    files.append(os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_manual_split_iterations_top3.svg'))
    files.append(os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_manual_split_iterations_top4.svg'))

    # Violin plot showing the rank distribution in the 10 iterations across different features on the x-axis
    #files.append(os.path.join(config['figures']['violin'], 'ranks_per_feature_10_iterations.svg'))

    # SHAP summary plot of all 10 iterations per sno_type (dot plot and bar plot)
    ###files.extend(expand(config['figures']['summary_shap_snotype'] + '{models3}_{sno_type}_test_set_100_background_scale_after_split_10_iterations.svg', **config))
    ###files.extend(expand(config['figures']['bar'] + '{models3}_{sno_type}_test_set_100_background_scale_after_split_10_iterations.svg', **config))

    # Logo of C box in expressed C/D box snoRNAs
    #files.append(os.path.join(config['figures']['logo'], 'expressed_c_prime_box.svg'))
    # Logo of H box in expressed H/ACA box snoRNAs
    #files.append(os.path.join(config['figures']['logo'], 'expressed_h_box.svg'))

    # Logo of C box in expressed C/D box snoRNAs with 3 nt flanking upstream and downstream
    #files.append(os.path.join(config['figures']['logo_w_flanking'], 'expressed_c_box_w_flanking_nt.svg'))
    # Logo of H box in expressed H/ACA box snoRNAs with 3 nt flanking upstream and downstream
    #files.append(os.path.join(config['figures']['logo_w_flanking'], 'expressed_h_box_w_flanking_nt.svg'))

    # Logo wo blank of C box in expressed C/D box snoRNAs
    files.append(os.path.join(config['figures']['logo'], 'expressed_c_prime_box_wo_blank.svg'))

    # Feature distribution
    '''
    files.extend(expand(config['figures']['density'] + '{numerical_features}.svg', **config))
    files.extend(expand(config['figures']['density'] + numerical, **config))
    #files.extend(expand(config['figures']['pairplot'] + cd_numerical, **config))
    #files.extend(expand(config['figures']['pairplot'] + haca_numerical, **config))
    files.extend(expand(config['figures']['bar'] + '{categorical_features}.svg', **config))
    files.extend(expand(config['figures']['bar_split_sno_type'] + '{categorical_features}_cd.svg', **config))
    files.extend(expand(config['figures']['bar_split_sno_type'] + '{categorical_features}_haca.svg', **config))
    files.extend(expand(config['figures']['density_split_sno_type'] + '{numerical_features}_cd.svg', **config))
    files.extend(expand(config['figures']['density_split_sno_type'] + '{numerical_features}_haca.svg', **config))
    #files.extend(expand(config['figures']['density_split_sno_type'] + '{numerical_features_scaled}_cd.svg', **config))
    #files.extend(expand(config['figures']['density_split_sno_type'] + '{numerical_features_scaled}_haca.svg', **config))
    files.extend(expand(config['figures']['bar_split_sno_type'] + '{intronic_features}_cd.svg', **config))
    files.extend(expand(config['figures']['bar_split_sno_type'] + '{intronic_features}_haca.svg', **config))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_sno_type.svg'))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_host_biotype.svg'))
    files.append(os.path.join(config['figures']['density_split_sno_type'], 'sno_mfe_length_normalized_haca.svg'))
    files.append(os.path.join(config['figures']['density_split_sno_type'], 'sno_mfe_length_normalized_cd.svg'))
    '''
    # Intron subgroup snoRNA feature distribution (small vs long intron)
    #files.extend(expand(config['figures']['density_split_sno_type'] + '{sno_type}_small_intron_{intron_group_feature}.svg', **config))
    #files.extend(expand(config['figures']['density_split_sno_type'] + '{sno_type}_long_intron_{intron_group_feature}.svg', **config))
    #files.append(os.path.join(config['figures']['donut'], 'abundance_status_intron_subgroup.svg'))



    # Mouse figures
    files.extend(expand(config['figures']['density'] + '{mouse_numerical_features}_abundance_status_mouse_{sno_type}.svg', **config))
    files.extend(expand(os.path.join(config['figures']['bar_split_sno_type'] + 'host_abundance_cutoff_mouse_{sno_type}.svg'), **config))
    files.append(os.path.join(config['figures']['violin'], 'models_accuracies_iterations_mouse.svg'))
    files.append(os.path.join(config['figures']['pie'], 'sno_number_per_confusion_value_mouse.svg'))
    files.extend(expand(os.path.join(config['figures']['pie'], 'sno_number_per_confusion_value_mouse_{models2}.svg'), **config))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_sno_type_mouse.svg'))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_host_biotype_mouse.svg'))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_sno_type_mouse_wo_dup.svg'))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_host_biotype_mouse_wo_dup.svg'))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_sno_type_mouse_no_dup.svg'))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_host_biotype_mouse_no_dup.svg'))
    files.append(os.path.join(config['figures']['violin'], "avg_tpm_per_confusion_value_mouse.svg"))
    files.extend(expand(os.path.join(config['figures']['violin'], "avg_tpm_per_confusion_value_mouse_{models2}.svg"), **config))
    files.extend(expand(os.path.join(config['figures']['pie'], 'sno_number_per_confusion_value_species_prediction_{models2}_{rs}.svg'), **config))
    #files.append(os.path.join(config['figures']['scatter'],
    #            'all_model_accuracies_cv_train_test_species_prediction.svg'))
    #files.append(os.path.join(config['figures']['scatter'],
    #            'all_model_accuracies_cv_train_test_species_prediction_rs.svg'))
    files.append(os.path.join(config['figures']['scatter'],
                'all_model_accuracies_cv_train_test_species_prediction_wo_dup_rs.svg'))
    files.append(os.path.join(config['figures']['scatter'],
                'all_model_accuracies_cv_train_test_species_prediction_no_dup_rs.svg'))
    #files.append(os.path.join(config['figures']['pie'], 'sno_number_per_confusion_value_species_prediction_log_reg_thresh.svg'))
    #files.append(os.path.join(config['figures']['scatter'],
    #                    'all_model_accuracies_cv_train_test_species_prediction_w_log_reg_thresh.svg'))
    files.append(os.path.join(config['figures']['donut'],
                            'confusion_value_host_biotype_mouse.svg'))
    files.extend(expand(os.path.join(config['figures']['donut'],
                        'confusion_value_host_biotype_mouse_{models2}_wo_dup.svg'), **config))
    files.extend(expand(os.path.join(config['figures']['donut'],
                        'confusion_value_host_biotype_mouse_{models2}_no_dup.svg'), **config))
    #files.append(os.path.join(config['figures']['violin'],
    #                "avg_tpm_per_confusion_value_mouse_log_reg_thresh.svg"))
    #files.append(os.path.join(config['figures']['scatter'],
    #                    'all_model_accuracies_cv_train_test_species_prediction_rs_w_log_reg_thresh.svg'))
    #files.append(os.path.join(config['figures']['scatter'],
    #                    'all_model_accuracies_cv_train_test_species_prediction_rs_w_log_reg_thresh_top3.svg'))
    files.append(os.path.join(config['figures']['scatter'],
                'all_model_accuracies_cv_train_test_species_prediction_wo_dup_top3_rs.svg'))
    files.append(os.path.join(config['figures']['donut'],
                        'confusion_value_host_biotype_mouse_top3_wo_dup_log_reg.svg'))  # this is not log_reg thresh


    # GTEx comparison vs TGIRT host abundance cutoff
    files.append(os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_manual_split_iterations_gtex_HG.svg'))
    files.append(os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_manual_split_gtex_HG.svg'))
    files.append(os.path.join(config['figures']['violin'], 'ranks_per_feature_manual_split_iterations_HG.svg'))
    files.append(os.path.join(config['figures']['venn'],
                        'gtex_tgirt_host_expressed_intersection.svg')),
    files.append(os.path.join(config['figures']['venn'],
                        'gtex_unpaired_tgirt_host_expressed_intersection.svg')),

    # Confusion value comparisons of predictive features
    #files.extend(expand(os.path.join(config['figures']['density_confusion_value'],
    #                        "{comparison_confusion_val}_{top_10_numerical_features}.svg"), **config))
    #files.extend(expand(os.path.join(config['figures']['density_confusion_value'],
    #                        "{comparison_confusion_val}_{top_10_numerical_features}_{sno_type}.svg"), **config))
    #files.extend(expand(os.path.join(config['figures']['bar_confusion_value'],
    #                       "{comparison_confusion_val}_{top_10_categorical_features}.svg"), **config))
    #files.extend(expand(os.path.join(config['figures']['bar_confusion_value'],
    #                        "{comparison_confusion_val}_{top_10_categorical_features}_{sno_type}.svg"), **config))

    # Comparison of TPMs for real FN vs real TP
#####    files.append(os.path.join(config['figures']['violin'], "FN_vs_TP_vs_TP_vs_FP_tpm.svg"))

    # Shap local value per confusion value across 10 iterations and per model
    #files.extend(expand(os.path.join(config['figures']['decision_plot_10_iterations'], '{confusion_value}_10_iterations_log_reg.svg'), **config))

    # FP vs TN in an expressed HG
    #files.extend(expand(os.path.join(config['figures']['density_confusion_value'],
    #                            "host_expressed_FP_vs_TN_{top_10_numerical_features}.svg"), **config))

    # Sashimi plot for C19orf48 (SNORD88B host gene)
    #files.extend(expand(os.path.join(config['figures']['sashimi'], '{tissue}_c19orf48.svg'), **config))

    # RBP par-clip analyses
    #files.append(os.path.join(config['figures']['density'], 'rbp_enrichment_combined.svg'))
    #files.append(os.path.join(config['figures']['density'], 'combined_rbp_score_confusion_value.svg'))

    # SHAP heatmaps
    #files.append(os.path.join(config['figures']['heatmap'], 'shap_all_models_iterations.svg'))
    #files.extend(expand(os.path.join(config['figures']['heatmap'], 'shap_{models2}_all_iterations.svg'), models2=config['models3']))
    #files.extend(expand(os.path.join(config['figures']['heatmap'], 'shap_{confusion_value}_all_models_iterations.svg'), **config))
    #files.extend(expand(os.path.join(config['figures']['heatmap'], 'shap_{confusion_value}_{models2}_all_iterations.svg'),
    #                    confusion_value=config['confusion_value'], models2=config['models3']))

    # SHAP top 3 features
    #files.extend(expand(os.path.join(config['figures']['bar_confusion_value'], 'shap_top_3_features_{confusion_value}.svg'), **config))

    # Multi HG different labels snoRNAs
#####    files.append(os.path.join(config['figures']['bar_confusion_value'], 'FP_vs_TN_multi_HG_different_labels.svg'))
    # Feature rank correlation across iterations
    #files.extend(expand(os.path.join(config['figures']['heatmap'],
    #            '{model_name}_feature_rank_correlation_10_iterations.svg'),
    #            model_name=['log_reg', 'svc', 'rf']))
    #files.extend(expand(os.path.join(config['figures']['heatmap'],
    #            '{model_name}_feature_rank_correlation_10_iterations_clustered.svg'),
    #            model_name=['log_reg', 'svc', 'rf']))
    #files.append(os.path.join(config['figures']['roc'], 'roc_curves_test_set_4_models.svg'))
    #files.append(os.path.join(config['figures']['roc'], 'roc_curves_test_set_5_models.svg'))
    #files.append(os.path.join(config['figures']['roc'], 'roc_curves_test_set_5_models_wo_clusters.svg'))
    #files.extend(expand(config['figures']['roc'] + 'roc_curves_test_set_5_models_wo_{feature_effect}.svg', **config))  # roc curve of performance without one or all top features
    #files.extend(expand(config['figures']['roc'] + 'roc_curves_test_set_5_models_only_{one_feature}.svg', **config))  # roc curve of performance with only one top predicting feature
    #files.append(os.path.join(config['figures']['roc'], 'roc_curves_test_set_5_models_scale_after_split.svg'))

    # Feature importance built-in sklearn
    #files.extend(expand(config['figures']['bar'] + '{models}_sklearn_feature_importance.svg', **config))
    #files.extend(expand(config['figures']['bar'] + '{models}_sklearn_feature_importance_wo_top_10_all.svg', **config))  #feature importance without the top_10_all features
    #files.extend(expand(config['figures']['bar'] + '{models}_sklearn_feature_importance_wo_top_11_all.svg', **config))

    # PCA plots
    #files.extend(expand(config['figures']['pca'] + "all_snoRNAs_{pca_hue}.svg", **config))
    #files.extend(expand(config['figures']['pca'] + "snoRNAs_test_{pca_hue}.svg", **config))

    # tSNE plots
    #files.extend(expand(config['figures']['t_sne'] + "all_snoRNAs_{pca_hue}.svg", **config))
    #files.extend(expand(config['figures']['t_sne'] + "snoRNAs_test_{pca_hue}.svg", **config))

    ## SHAP summary plots for all snoRNAs in test set or per sno_type or hg_biotype
    #files.extend(expand(config['figures']['summary_shap'] + '{models}_all_features_test_set.svg', **config)) # "Train" the explainer on all the training set and "explain" on test set
    #files.extend(expand(config['figures']['summary_shap'] + '{models}_all_features_test_set_100_background.svg', **config)) # "Train" the explainer on 100 examples of the training set and "explain" on test set
    #files.extend(expand(config['figures']['summary_shap_snotype'] + '{models}_{sno_type}_all_features_test_set_100_background.svg', ** config))  # Same, but per sno_type
    #files.extend(expand(config['figures']['summary_shap_hg_biotype'] + '{models}_{hg_biotype}_all_features_test_set_100_background.svg', **config))  # Same, but per intronic vs intergenic snoRNAs
    ##files.extend(expand(config['figures']['decision_plot'] + 'ENSG00000226572_{models}_all_features_test_set_100_background.svg', **config))  # Decision plot for the specific ENSG00000226572 snoRNA
    ##files.extend(expand(config['figures']['decision_plot'] + 'ENSG00000221116_{models}_all_features_test_set_100_background.svg', **config))  # Decision plot for the specific ENSG00000221116 snoRNA
    ##files.extend(expand(config['figures']['decision_plot'] + 'ENSG00000212498_{models}_all_features_test_set_100_background.svg', **config))  # Decision plot for the specific ENSG00000212498 snoRNA

    ## SHAP summary plots for all snoRNAs in test set per sno_type or hg_biotype (this time with scale of values after split of datasets)
    #files.extend(expand(config['figures']['summary_shap_snotype'] + '{models2}_{sno_type}_test_set_100_background_scale_after_split.svg', ** config))
    #files.extend(expand(config['figures']['summary_shap_hg_biotype'] + '{models2}_{hg_biotype}_test_set_100_background_scale_after_split.svg', **config))

    ## Upset plots
    #files.extend(expand(config['figures']['upset'] + '{confusion_value}_all_models.svg', **config))
    #files.extend(expand(config['figures']['upset'] + '{confusion_value}_all_models_scale_after_split.svg', **config))
    #files.append(os.path.join(config['figures']['upset'], 'top_5_features_intersection_all_models_scale_after_split.svg'))

    ## Donut charts linked to upset plot top 5 features
    #files.append(os.path.join(config['figures']['donut'], 'top_feature_intersection_rank_percent.svg'))

    ## Connected scatter plot of accuracies of all models in CV, train and test sets
    #files.append(os.path.join(config['figures']['scatter'], 'all_model_accuracies_cv_train_test.svg'))

    ## Violin plot of the predictive rank of all feature (distribution across all 4 models)
    #files.append(os.path.join(config['figures']['violin'], 'ranks_per_feature.svg'))

    ## Pair plot per sno_type of the predictive rank of all feature (distribution across all 4 models)
    #files.append(os.path.join(config['figures']['pairplot'], 'top_5_features_cd.svg'))
    #files.append(os.path.join(config['figures']['pairplot'], 'top_5_features_haca.svg'))

    ## Feature importance global shap bar plot
    #files.extend(expand(config['figures']['bar'] + '{models}_global_feature_importance.svg', **config))
    #files.extend(expand(config['figures']['bar'] + '{models2}_global_feature_importance_scale_after_split.svg', **config))

    ## Heatmap shap values across test set
    #files.extend(expand(config['figures']['heatmap'] + '{models2}_shap_heatmap.svg', **config))

    ## Forgi structure plot for SNORA77B
    #files.append(os.path.join(config['figures']['forgi_structure'],  'SNORA77B_terminal_stem.svg'))

    ## Bivariate density plot showing conservation and structure stability of H/ACA snoRNAs
    #files.append(os.path.join(config['figures']['bivariate_density'], "HACA_abundance_status_conservation_mfe.svg"))

    ## Violin plot to compare the abundance of expressed C/D and H/ACA snoRNAs
    #files.append(os.path.join(config['figures']['violin'], 'tpm_comparison_sno_type.svg'))

    ## Violin plot to compare the abundance of expressed H/ACA and C/D with or without terminal stem
    #files.append(os.path.join(config['figures']['violin'], 'cd_stem_comparison_tpm.svg'))
    #files.append(os.path.join(config['figures']['violin'], 'haca_stem_comparison_tpm.svg'))

    return files
