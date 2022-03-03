import os
include: "figures_model_output.smk"

""" Cluster snoRNAs based on their SHAP values to see if any subgroup immerges
    -per model or per iteration?
    -per label (expressed vs not expressed)
    -per confusion value (TP, TN, FP, FN)"""

rule heatmap_shap_all_models_iterations:
    """ Cluster snoRNAs based on the shap values of their features across all
        iterations and models (log_reg, svc and rf)"""
    input:
        shap_values = expand(rules.get_all_shap_values.output.shap, iteration=config['iteration'], models2=config['models3'])
    output:
        heatmap = os.path.join(config['figures']['heatmap'], 'shap_all_models_iterations.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/heatmap_shap_all_models_iterations.py"

rule heatmap_shap_per_model_all_iterations:
    """ Cluster snoRNAs based on the shap values of their features across all
        iterations per model (log_reg, svc or rf)"""
    input:
        shap_values = expand(rules.get_all_shap_values.output.shap, iteration=config['iteration'], models2=config['models3'])
    output:
        heatmap = expand(os.path.join(config['figures']['heatmap'], 'shap_{models2}_all_iterations.svg'), models2=config['models3'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/heatmap_shap_per_model_all_iterations.py"

rule heatmap_shap_per_confusion_value_all_models_iterations:
    """ Cluster snoRNAs based on the shap values of their features across all
        iterations/models per confusion_value"""
    input:
        shap_values = expand(rules.get_all_shap_values.output.shap, iteration=config['iteration'], models2=config['models3']),
        confusion_value_df = rules.real_confusion_value_df.output.real_confusion_value_df
    output:
        heatmap = os.path.join(config['figures']['heatmap'], 'shap_{confusion_value}_all_models_iterations.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/heatmap_shap_per_confusion_value_all_models_iterations.py"

rule heatmap_shap_per_confusion_value_per_model_all_iterations:
    """ Cluster snoRNAs based on the shap values of their features across all
        iterations/models per confusion_value"""
    input:
        shap_values = expand(rules.get_all_shap_values.output.shap, iteration=config['iteration'], allow_missing=True),
        confusion_value_df = rules.real_confusion_value_df.output.real_confusion_value_df
    output:
        heatmap = os.path.join(config['figures']['heatmap'], 'shap_{confusion_value}_{models2}_all_iterations.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/heatmap_shap_per_confusion_value_per_model_all_iterations.py"

rule bar_shap_top_3_features_per_confusion_value:
    """ Find the relative number of times per confusion value that a feature is the
        top 1, 2 or 3 based on the highest absolute value of SHAP across models
        and iterations. Drop duplicate snoRNAs if they have the same top1/2/3
        and iteration. Return a bar chart (and the corresponding df) for each
        confusion value."""
    input:
        shap_values = expand(rules.get_all_shap_values.output.shap, iteration=config['iteration'], models2=config['models3']),
        confusion_value_df = rules.real_confusion_value_df.output.real_confusion_value_df
    output:
        bar = os.path.join(config['figures']['bar_confusion_value'], 'shap_top_3_features_{confusion_value}.svg'),
        df = os.path.join(config['path']['real_confusion_value'], 'shap_top_3_features_{confusion_value}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/bar_shap_top_3_features_per_confusion_value.py"
