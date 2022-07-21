import os
include: "figures_model_output.smk"
include: "cv_train_test_manual_split.smk"
include: "downloads.smk"

rule filter_bam:
    """ Filter bam to only keep the reads in th region of the SNORD88B host gene
        C19orf48. """
    input:
        bams = os.path.join(config['path']['bams'], '{tissue}/Aligned.sortedByCoord.out.bam')
    output:
        bam_index = os.path.join(config['path']['bams'], '{tissue}/Aligned.sortedByCoord.out.bam.bai'),
        filtered_bams = os.path.join(config['path']['filtered_bams'], '{tissue}_c19orf48.bam'),
        bam_index_filtered = os.path.join(config['path']['filtered_bams'], '{tissue}_c19orf48.bam.bai')
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input.bams} {output.bam_index} && "
        "samtools view {input.bams} '19:50797704-50804929' -b > {output.filtered_bams} && "
        "samtools index {output.filtered_bams} {output.bam_index_filtered}"

rule sashimi_snord88_host_gene:
    """ Create a sashimi plot to show the splicing patterns of C19orf48
        (SNORD88B host gene) in each tissue."""
    input:
        bams = rules.filter_bam.output.filtered_bams,
        sashimi_script = rules.sashimi_script_download.output.sashimi_script,
        gtf = config['path']['gtf']
    output:
        sashimi = os.path.join(config['figures']['sashimi'], '{tissue}_c19orf48.svg')
    conda:
        "../envs/sashimi.yaml"
    shell:
        "./{input.sashimi_script} -b {input.bams} -c '19:50797704-50804929' -o {output.sashimi} -g {input.gtf} -F 'svg'"

rule multi_HG_different_label_snoRNAs:
    """ Find all the host genes (HG) containing multiple snoRNAs and that do not
        have the same label (i.e. all HG containing at least 1 expressed snoRNA
        and 1 not expressed snoRNA)."""
    input:
        host_df = config['path']['host_gene_df'],
        feature_df = rules.merge_feature_df.output.feature_df
    output:
        multi_HG_different_label_snoRNAs = config['path']['multi_HG_different_label_snoRNAs']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/multi_HG_different_label_snoRNAs.py"

rule bar_FP_vs_TN_multi_HG_different_labels:
    """ Create a bar chart showing the proportion of snoRNAs included within a
        HG that has snoRNAs with different labels for FP vs TN."""
    input:
        confusion_value_df = expand(rules.real_confusion_value_df.output.real_confusion_value_df, confusion_value=['FP', 'TN']),
        multi_HG_df = rules.multi_HG_different_label_snoRNAs.output.multi_HG_different_label_snoRNAs
    output:
        bar = os.path.join(config['figures']['bar_confusion_value'], 'FP_vs_TN_multi_HG_different_labels.svg')
    params:
        color_dict = config['colors_complex']['multi_HG_labels']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/bar_FP_vs_TN_multi_HG_different_labels.py"

rule decision_plot_interesting_snoRNAs:
    """ Create a decision plot for intersting snoRNAs (ex: SNORA77B, SNORD86) (for all models)."""
    input:
        X_test = expand(rules.fill_na_feature_scaling_after_manual_split.output.test, manual_iteration=config['manual_iteration']),
        shap_values = expand(rules.get_all_shap_values_manual_split.output.shap, manual_iteration=config['manual_iteration'], allow_missing=True),
        expected_values = expand(rules.get_all_shap_values_manual_split.output.expected_value, manual_iteration=config['manual_iteration'], allow_missing=True)
    output:
        decision_plot = os.path.join(config['figures']['decision_plot_interesting_snoRNAs'], '{interesting_sno_ids}_{models2}.svg')
    params:
        colors_dict = config['colors_complex']['abundance_cutoff_2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/decision_plot_interesting_snoRNAs.py"

rule decision_plot_multi_HG_snoRNAs:
    """ Create a decision plot containing all snoRNAs in GAS5 and SNHG17 (multi HG with
        snoRNAs with different labels)."""
    input:
        shap_values = expand(rules.get_all_shap_values_manual_split.output.shap,
                        manual_iteration=config['manual_iteration'], allow_missing=True),
        expected_values = expand(rules.get_all_shap_values_manual_split.output.expected_value,
                        manual_iteration=config['manual_iteration'], allow_missing=True)
    output:
        decision_plot = os.path.join(config['figures']['decision_plot'], '{multi_HG_diff_label}_snoRNAs_{models2}.svg')
    params:
        sno_ids = lambda wildcards: config['multi_HG_sno_ids'][wildcards.multi_HG_diff_label],
        colors_dict = config['colors_complex']['abundance_cutoff_2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/decision_plot_multi_HG_snoRNAs.py"

rule pie_donut_multi_HG_snoRNAs:
    """ Create different pie chart showing mono vs multi_HG, multi-HG with same
        vs different snoRNA labels and a donut chart showing multi_hg with
        different labels that are 50/50 expressed-not_expressed vs debalanced
        (outer donut) and the confusion_value proportions (inner_value)."""
    input:
        df = rules.merge_feature_df.output.feature_df,
        host_df = config['path']['host_gene_df'],
        multi_HG_different_label_snoRNAs = rules.multi_HG_different_label_snoRNAs.output.multi_HG_different_label_snoRNAs,
        fake_dependency = expand(rules.confusion_matrix_f1_scale_after_manual_split.output.info_df,
                                manual_iteration=config['manual_iteration'], models2='log_reg', allow_missing=True),
        fake_dependency2 = expand(rules.upset_models_confusion_scale_after_manual_split.output.upset, **config) 
    output:
        pie = os.path.join(config['figures']['pie'], 'mono_vs_multi_HG.svg'),
        donut = os.path.join(config['figures']['donut'], 'multi_HG_diff_sno_labels_proportion_and_confusion_value.svg')
    params:
        merged_confusion_values = expand("results/tables/confusion_matrix_f1/merged_confusion_matrix_{manual_iteration}.tsv",
                                manual_iteration=config['manual_iteration']),
        mono_vs_multi_HG_colors = config['colors_complex']['mono_vs_multi_HG'],
        multi_HG_labels_colors = config['colors_complex']['multi_HG_labels'],
        multi_HG_same_labels_proportion_colors = config['colors_complex']['multi_HG_same_labels_proportion'],
        multi_HG_diff_labels_proportion_colors = config['colors_complex']['multi_HG_diff_labels_proportion']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/pie_donut_multi_HG_snoRNAs.py"

rule hbar_nb_sno_per_confusion_val_HG:
    """ Create a stacked horizontal bar chart where each bar corresponds to the
        stacked number of snoRNAs in each confusion value. Sort the bars
        according to if the snoRNAs have the same or different labels among the
        HG and according to the expression categories (see pie_donut_multi_HG_snoRNAs)"""
    input:
        sno_per_confusion_value = rules.regroup_sno_confusion_value_manual_split.output.sno_per_confusion_value,
        host_df = config['path']['host_gene_df'],
        multi_HG_different_label_snoRNAs = rules.multi_HG_different_label_snoRNAs.output.multi_HG_different_label_snoRNAs,
        feature_df = rules.merge_feature_df.output.feature_df
    output:
        hbar = os.path.join(config['figures']['hbar'], 'nb_sno_per_confusion_val_multi_HG.svg')
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/hbar_nb_sno_per_confusion_val_HG.py"
