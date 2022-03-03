import os
import subprocess as sp

include: "downloads.smk"
include: "data_processing.smk"
include: "figures_model_output.smk"

""" Deal with PAR-CLIP datasets."""


rule merge_bed_peaks:
    """ For each bed obtained from the PAR-CLIP datasets, merge peaks on same
        strand that overlap by summing their enrichment value."""
    input:
        input_beds = rules.par_clip_download.output
    output:
        merge_beds = expand(os.path.join(config['path']['par_clip'], '{rbp}_merged.bed'),
                            rbp=['NOP58_repA', 'NOP58_repB', 'NOP56', 'FBL', 'FBL_mnase', 'DKC1'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/merge_bed_peaks.py"


rule format_snoRNA_bed:
    """ Remove last 6 columns of the snoRNA bed file."""
    input:
        snoRNA_bed = rules.gtf_to_bed.output.all_sno_bed
    output:
        formated_bed = config['path']['all_sno_bed_formated']
    shell:
        "cut -f1-6 {input.snoRNA_bed} > {output.formated_bed}"

rule sort_par_clip_bed:
    """ Sort alphabetically PAR-CLIP beds."""
    input:
        beds = rules.merge_bed_peaks.output.merge_beds
    output:
        sorted_beds = expand(os.path.join(config['path']['par_clip'], '{rbp}_merged_sorted.bed'),
                            rbp=['NOP58_repA', 'NOP58_repB', 'NOP56', 'FBL', 'FBL_mnase', 'DKC1'])
    run:
        for i, bed in enumerate(input.beds):
            output_i = output.sorted_beds[i]
            sp.call("""awk -v OFS="\t" '{print $1,$2,$3,$4,".",$5}' """+bed+""" | sort -k1,1 -k2n > """+output_i, shell=True)

rule map_bed_peaks_to_sno:
    """ For each merged bed obtained from merge_bed_peaks, map the enrichment
        value of overlapping peaks to each snoRNA (return the sum of all peaks
        overlapping a snoRNA per RBP)."""
    input:
        merged_beds = rules.sort_par_clip_bed.output.sorted_beds,
        snoRNA_bed = rules.format_snoRNA_bed.output.formated_bed
    output:
        mapped_snoRNA_bed = expand(os.path.join(config['path']['par_clip'], '{rbp}_mapped_to_snoRNAs.bed'),
                            rbp=['NOP58_repA', 'NOP58_repB', 'NOP56', 'FBL', 'FBL_mnase', 'DKC1'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/map_bed_peaks_to_sno.py"

rule rbp_enrichment_density:
    input:
        all_feature_df = config['path']['feature_df'],
        beds = rules.map_bed_peaks_to_sno.output.mapped_snoRNA_bed
    output:
        density = expand(os.path.join(config['figures']['density'], '{rbp}_enrichment.svg'), rbp=['NOP58_repA', 'NOP58_repB', 'NOP56', 'FBL', 'FBL_mnase', 'DKC1']),
        density_combined = os.path.join(config['figures']['density'], 'rbp_enrichment_combined.svg'),
        combined_rbp_score_df = config['path']['combined_rbp_score_df']
    params:
        color_dict = config['colors_complex']['label']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/rbp_enrichment_density.py"

rule rbp_enrichment_density_confusion_value:
    input:
        confusion_value_dfs = expand(rules.real_confusion_value_df.output.real_confusion_value_df, **config),
        combined_rbp_score_df = rules.rbp_enrichment_density.output.combined_rbp_score_df
    output:
        density = os.path.join(config['figures']['density'], 'combined_rbp_score_confusion_value.svg'),
        dfs = expand(os.path.join(config['path']['real_confusion_value'], '{confusion_value}_w_rbp_features.tsv'), **config)
    params:
        color_dict = config['colors_complex']['confusion_value']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/rbp_enrichment_density_confusion_value.py"
