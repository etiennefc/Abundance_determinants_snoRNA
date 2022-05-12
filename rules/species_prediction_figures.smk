import os

rule density_features_species:
    """ Generate a density plot of the distribution for all numerical input
        features. The separation is based on the predicted abundance status. """
    input:
        df = rules.predict_species_snoRNA_label.output.predicted_label_df,
        snoRNA_type_df = rules.find_species_snoRNA_type.output.snoRNA_type_df
    output:
        density_features = os.path.join(config['figures']['density'],
                            '{species_numerical_features}_abundance_status_{species}_{sno_type}.svg')
    params:
        hue_color = config['colors_complex']['label']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/density_features_species.py"

rule bar_host_expressed_species:
    """ Generate a bar chart of all predicted expressed vs non_expressed snoRNAs
        according to the host_abundance_cutoff feature (separately for C/D and
        H/ACA snoRNAs)."""
    input:
        df = rules.predict_species_snoRNA_label.output.predicted_label_df,
        snoRNA_type_df = rules.find_species_snoRNA_type.output.snoRNA_type_df
    output:
        bar = os.path.join(config['figures']['bar_split_sno_type'],
                            'host_abundance_cutoff_{species}_{sno_type}.svg')
    params:
        hue_color = config['colors_complex']['abundance_cutoff_host']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/bar_categorical_species.py"

rule donut_predicted_label_sno_type_species:
    """ Generate a donut chart of the number and % of predicted expressed vs
        not expressed snoRNAs (outer donut) and per sno_type (inner donut) for
        species snoRNAs."""
    input:
        df = rules.predict_species_snoRNA_label.output.predicted_label_df,
        snoRNA_type_df = rules.find_species_snoRNA_type.output.snoRNA_type_df
    output:
        donut = os.path.join(config['figures']['donut'],
                            'predicted_abundance_status_sno_type_{species}.svg')
    params:
        label_colors = config['colors_complex']['abundance_cutoff_2'],
        sno_type_colors = config['colors_complex']['sno_type']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/donut_labels_sno_type_species.py"

rule donut_label_host_biotype_species:
    """ Generate a donut chart of the number and % of predicted expressed vs
        not expressed snoRNAs (outer donut) and per host biotype (inner donut)
        for species snoRNAs."""
    input:
        df = rules.predict_species_snoRNA_label.output.predicted_label_df,
        host_biotype_df = rules.find_species_snoRNA_HG.output.species_snoRNA_HG
    output:
        donut = os.path.join(config['figures']['donut'],
                            'predicted_abundance_status_host_biotype_{species}.svg')
    params:
        label_colors = config['colors_complex']['abundance_cutoff_2'],
        host_biotype_colors = config['colors_complex']['host_biotype2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/donut_labels_host_biotype_species.py"

rule bar_ab_status_prediction_species:
    """ Generate a stacked bar chart of all predicted expressed vs not_expressed
        snoRNAs per species and also actual abundance status proportion for
        human and mouse. Add also the total number of snoRNAs on top of each
        bar."""
    input:
        dfs = expand(rules.predict_species_snoRNA_label.output.predicted_label_df, species=species),
        mouse_labels = rules.merge_features_label_mouse.output.feature_df,
        human_labels = rules.merge_feature_df.output.feature_df
    output:
        bar = os.path.join(config['figures']['bar'],
                            'ab_status_prediction_all_species.svg')
    params:
        hue_color = config['colors_complex']['abundance_cutoff_2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/bar_ab_status_prediction_species.py"
