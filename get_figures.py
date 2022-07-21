#!/usr/bin/python3
from snakemake.io import expand
import os
def get_figures_path(config):
    """Return list of figures (and their path) to generate from config"""
    files = []
    numerical = "{numerical_features}_{feature_hue}.svg"
    cd_numerical = "numerical_features_cd.svg"
    haca_numerical = "numerical_features_haca.svg"

    # ROC curve and connected scatter plot of average accuracies of all models in CV, train and test sets across 10 manual split iterations
    files.append(os.path.join(config['figures']['roc'],
                    'roc_curves_test_set_5_models_scale_after_manual_split.svg'))
    files.append(os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_manual_split_iterations.svg'))

    # Violin plot of feature rank across 10 manual split iterations and upset plot of confusion value
    files.append(os.path.join(config['figures']['violin'], 'ranks_per_feature_manual_split_iterations.svg'))
    files.extend(expand(os.path.join(config['figures']['upset'],
                    '{confusion_value}_all_models_scale_after_manual_split.svg'), **config))

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

    # Bar plots to compare per confusion_value the top 10 categorical features distributions
    files.extend(expand(os.path.join(config['figures']['bar_confusion_value'],
                                "confusion_value_comparison_{top_10_categorical_features}.svg"), **config))

    # Violin plot to show the log2(avg TPM) per confusion value
    files.append(os.path.join(config['figures']['violin'], "avg_tpm_per_confusion_value.svg"))

    # Horizontal stacked bar chart to show the number of sno per HG (per confusion value)
    files.append(os.path.join(config['figures']['hbar'], 'nb_sno_per_confusion_val_multi_HG.svg'))

    # ROC curve of models trained with top1 (only hamming box score), top3 or top4 features across 10 manual iterations
    files.append(os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_manual_split_top3.svg'))
    files.append(os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_manual_split_top4.svg'))
    files.append(os.path.join(config['figures']['roc'],
                        'roc_curves_test_set_5_models_scale_after_split_10_iterations_only_hamming.svg'))

    # Connected scatter plot of average accuracies of all models in CV, train and test sets across 10 manual iterations of models trained on only hamming box score, top3 or top4 features
    files.append(os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_manual_split_iterations_top3.svg'))
    files.append(os.path.join(config['figures']['scatter'],
                        'all_model_accuracies_cv_train_test_manual_split_iterations_top4.svg'))
    files.append(os.path.join(config['figures']['scatter'], 
                        'all_model_accuracies_cv_train_test_10_iterations_only_hamming.svg'))

    # Logo wo blank of C box in expressed C/D box snoRNAs
    files.append(os.path.join(config['figures']['logo'], 'expressed_c_prime_box_wo_blank.svg'))

    # Feature distribution
    files.extend(expand(config['figures']['bar'] + '{categorical_features}.svg', **config))
    files.extend(expand(config['figures']['bar_split_sno_type'] + '{categorical_features}_cd.svg', **config))
    files.extend(expand(config['figures']['bar_split_sno_type'] + '{categorical_features}_haca.svg', **config))
    files.extend(expand(config['figures']['density_split_sno_type'] + '{numerical_features}_cd.svg', **config))
    files.extend(expand(config['figures']['density_split_sno_type'] + '{numerical_features}_haca.svg', **config))
    files.extend(expand(config['figures']['bar_split_sno_type'] + '{intronic_features}_cd.svg', **config))
    files.extend(expand(config['figures']['bar_split_sno_type'] + '{intronic_features}_haca.svg', **config))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_sno_type.svg'))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_host_biotype.svg'))
    
    # Mouse figures
    files.extend(expand(config['figures']['density'] + '{mouse_numerical_features}_abundance_status_mouse_{sno_type}.svg', **config))
    files.extend(expand(os.path.join(config['figures']['bar_split_sno_type'] + 'host_abundance_cutoff_mouse_{sno_type}.svg'), **config))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_sno_type_mouse_wo_dup.svg'))
    files.append(os.path.join(config['figures']['donut'], 'abundance_status_host_biotype_mouse_wo_dup.svg'))
    files.append(os.path.join(config['figures']['scatter'],
                'all_model_accuracies_cv_train_test_species_prediction_wo_dup_rs.svg'))
    files.extend(expand(os.path.join(config['figures']['donut'],
                        'confusion_value_host_biotype_mouse_{models2}_wo_dup.svg'), **config))

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


    return files
