import os
include: "cv_train_test_manual_split_gtex_HG.smk"
include: "feature_normalization.smk"
include: "gtex_HG_cutoff.smk"

rule create_local_env:
    """ Import matplotlib and seaborn in local snakemake environment. Note that the log
        will be empty, it serves only to link this rule to a decoy output"""
    output:
        log_create_env = "log/create_local_env.log"
    shell:
        """
        conda install -c conda-forge --force-reinstall seaborn
        conda install -c conda-forge --force-reinstall matplotlib
        pip install pca
        &> {output}
        """

rule pie_labels:
    """ Generate a pie chart of the number and % of expressed vs not expressed
        snoRNAs for all annotated snoRNAs."""
    input:
        df = config['path']['feature_df']
    output:
        pie = os.path.join(config['figures']['pie'], 'abundance_status.svg')
    params:
        colors = config['colors_complex']['abundance_cutoff_2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/pie.py"

rule donut_labels_sno_type:
    """ Generate a donut chart of the number and % of expressed vs not expressed
        snoRNAs (outer donut) and per sno_type (inner donut)."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        donut = os.path.join(config['figures']['donut'],
                            'abundance_status_sno_type.svg')
    params:
        label_colors = config['colors_complex']['abundance_cutoff_2'],
        sno_type_colors = config['colors_complex']['sno_type']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/donut_labels_sno_type.py"

rule donut_labels_host_biotype:
    """ Generate a donut chart of the number and % of expressed vs not expressed
        snoRNAs (outer donut) and per host biotype (inner donut)."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        donut = os.path.join(config['figures']['donut'],
                            'abundance_status_host_biotype.svg')
    params:
        label_colors = config['colors_complex']['abundance_cutoff_2'],
        host_biotype_colors = config['colors_complex']['host_biotype2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/donut_labels_host_biotype.py"

rule density_features_simple:
    """ Generate a simple density plot of the distribution for all numerical
        input features without any separation. """
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        density_features_simple = os.path.join(config['figures']['density'],
                                    '{numerical_features}.svg')
    params:
        simple_color = config['colors_complex']['simple']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/density_features_simple.py"

rule density_features:
    """ Generate a density plot of the distribution for all numerical input
        features. The separation can be done on either the sno_type,
        abundance_cutoff or abundance_cutoff_2. """
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        density_features = os.path.join(config['figures']['density'],
                            '{numerical_features}_{feature_hue}.svg')
    params:
        hue_color = lambda wildcards: config['colors_complex'][wildcards.feature_hue]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/density_features.py"

rule density_features_split:
    """ Generate a density plot of the distribution for all numerical input
        features by separating by sno_type and abundance_cutoff_2 (so 1 density
        plot for C/D and 1 for H/ACA per feature). """
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        density_features_cd = os.path.join(config['figures']['density_split_sno_type'],
                            '{numerical_features}_cd.svg'),
        density_features_haca = os.path.join(config['figures']['density_split_sno_type'],
                            '{numerical_features}_haca.svg')
    params:
        hue_color = config['colors_complex']['abundance_cutoff_2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/density_features_split.py"

rule density_normalized_mfe_split:
    """ Generate a density plot of the distribution of the sno mfe normalized
        by sno length by separating by sno_type and abundance_cutoff_2 (so 1
        density plot for C/D and 1 for H/ACA per feature)."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        density_features_cd = os.path.join(config['figures']['density_split_sno_type'],
                            'sno_mfe_length_normalized_cd.svg'),
        density_features_haca = os.path.join(config['figures']['density_split_sno_type'],
                            'sno_mfe_length_normalized_haca.svg')
    params:
        hue_color = config['colors_complex']['abundance_cutoff_2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/density_normalized_mfe_split.py"

#rule density_scaled_features_split:
#    """ Generate a density plot of the distribution for all SCALED numerical
#        features by separating by sno_type and abundance_cutoff_2 (so 1 density
#        plot for C/D and 1 for H/ACA per feature). """
#    input:
#        df = ****scaled_feature_df scale after split, but which iteration of cv-train-test split?****
#    output:
#        density_features_cd = os.path.join(config['figures']['density_split_sno_type'],
#                            '{numerical_features_scaled}_cd.svg'),
#        density_features_haca = os.path.join(config['figures']['density_split_sno_type'],
#                            '{numerical_features_scaled}_haca.svg')
#    params:
#        hue_color = config['colors_complex']['abundance_cutoff_2']
#    conda:
#        "../envs/python.yaml"
#    script:
#        "../scripts/python/graphs/density_scaled_features_split.py"

rule density_intron_groups_sno_type_features:
    """ Split intronic snoRNAs into small (<5000 nt) vs long (>=5000 nt) introns
        and create for each sno_type and each of these intron subgroup a density
        plot. The density plot shows the difference between expressed and not
        expressed snoRNAs across a variety of numerical features."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        density_small = os.path.join(config['figures']['density_split_sno_type'],
                            '{sno_type}_small_intron_{intron_group_feature}.svg'),
        density_long = os.path.join(config['figures']['density_split_sno_type'],
                            '{sno_type}_long_intron_{intron_group_feature}.svg')
    params:
        hue_color = config['colors_complex']['abundance_cutoff_2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/density_intron_groups_sno_type_features.py"

rule donut_labels_intron_subgroup:
    """ Create a donut chart of the expressed and not expressed snoRNAs (inner
        donut) per intron length subgroup."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        donut = os.path.join(config['figures']['donut'],
                            'abundance_status_intron_subgroup.svg')
    params:
        label_colors = config['colors_complex']['abundance_cutoff_2'],
        intron_subgroup_colors = config['colors_complex']['intron_subgroup']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/donut_labels_intron_subgroup.py"


rule pairplot_features:
    """ Generate a pair plot of all numerical input features with a hue of
        abundance_cutoff_2. Generate also the same pair plot, but for either C/D
        or H/ACA snoRNAs separately."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        pairplot = os.path.join(config['figures']['pairplot'],
                            'numerical_features.svg'),
        pairplot_cd = os.path.join(config['figures']['pairplot'],
                            'numerical_features_cd.svg'),
        pairplot_haca = os.path.join(config['figures']['pairplot'],
                            'numerical_features_haca.svg')
    params:
        hue_color = config['colors_complex']['abundance_cutoff_2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/pairplot_numerical.py"

rule bar_categorical_features:
    """ Generate a bar chart of all expressed vs non_expressed snoRNAs
        (abundance_cutoff_2) for the categorical features. Generate also the
        same bar graph but separately for C/D and H/ACA snoRNAs."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        bar_categorical_features = os.path.join(config['figures']['bar'],
                            '{categorical_features}.svg'),
        bar_categorical_features_cd = os.path.join(config['figures']['bar_split_sno_type'],
                            '{categorical_features}_cd.svg'),
        bar_categorical_features_haca = os.path.join(config['figures']['bar_split_sno_type'],
                            '{categorical_features}_haca.svg')
    params:
        hue_color = lambda wildcards: config['colors_complex'][wildcards.categorical_features]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/bar_categorical.py"

rule bar_large_scale_features:
    """ Generate a grouped bar chart for ranges of values for numerical features
        that have large data scale called intronic_features (e.g. dist_to_bp)
        with a hue of expressed vs non_expressed snoRNAs (abundance_cutoff_2)
        (and one graph per snoRNA type)."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        bar_cd = os.path.join(config['figures']['bar_split_sno_type'],
                        '{intronic_features}_cd.svg'),
        bar_haca = os.path.join(config['figures']['bar_split_sno_type'],
                            '{intronic_features}_haca.svg')
    params:
        hue_color = config['colors_complex']['abundance_cutoff_2']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/bar_intronic_features.py"

rule venn_host_abundance_cutoff_tgirt_gtex:
    """ Generate a Venn diagram of all intronic snoRNAs based on the feature
        host_expressed defined by TGIRT-Seq datasets or by GTEx datasets (one
        for host_expressed and one for host_not_expressed)."""
    input:
        gtex_df = rules.one_hot_encode_before_split_gtex_HG.output.one_hot_encoded_df,
        tgirt_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df
    output:
        venn_host_expressed = os.path.join(config['figures']['venn'],
                        'gtex_tgirt_host_expressed_intersection.svg'),
        venn_host_not_expressed = os.path.join(config['figures']['venn'],
                            'gtex_tgirt_host_not_expressed_intersection.svg')
    conda:
        "../envs/venn_diagram.yaml"
    script:
        "../scripts/python/graphs/venn_host_abundance_cutoff_tgirt_gtex.py"

rule venn_host_abundance_cutoff_tgirt_gtex_unpaired:
    """ Generate a Venn diagram of all intronic snoRNAs based on the feature
        host_expressed defined by TGIRT-Seq datasets or by unpaied tissue GTEx
        datasets (one for host_expressed and one for host_not_expressed)."""
    input:
        gtex_df = rules.one_hot_encode_before_split_gtex_HG_unpaired.output.one_hot_encoded_df,
        tgirt_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df
    output:
        venn_host_expressed = os.path.join(config['figures']['venn'],
                        'gtex_unpaired_tgirt_host_expressed_intersection.svg'),
        venn_host_not_expressed = os.path.join(config['figures']['venn'],
                            'gtex_unpaired_tgirt_host_not_expressed_intersection.svg')
    conda:
        "../envs/venn_diagram.yaml"
    script:
        "../scripts/python/graphs/venn_host_abundance_cutoff_tgirt_gtex.py"
