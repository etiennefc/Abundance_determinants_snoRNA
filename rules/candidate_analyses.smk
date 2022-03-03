import os
include: "figures_model_output.smk"
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

rule decision_plot_snora77b:
    """ Create a decision plot for all models and iterataions that SNORA77b is
        present in."""
    input:
        shap_val = expand(rules.get_all_shap_values.output.shap, iteration=config['iteration'], models2=config['models3'])
    output:
        fake_log = 'log/fake_log_snora77b.log'
    params:
        decision_plot_dir = config['figures']['decision_plot_snora77b']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/decision_plot_snora77b.py"
