import os

include: "cv_train_test_manual_split_top4.smk"

""" Predict SNORA81s abundance status for Laurence's analyses where SNORA81 is
    either WT or mutated in different regions of the snoRNA. The feature_df was
    generated in local (see rules/snora81_overexpression.smk)."""


rule predict_snora81s_label:
    """ Predict the abundance status of all SNORA81 (WT or mutants) based on the
        top4 features (combined_box_hamming. sno_mfe, terminal_stem_mfe and
        host_expressed). We scale the feature_df using the same parameters (mean
        and stdev) used to scale each training set."""
    input:
        feature_df = 'results/tables/all_features_snora81_top4.tsv',
        human_snoRNA_feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        pickled_trained_model = rules.train_models_scale_after_manual_split_top4.output.pickled_trained_model
    output:
        predicted_label_df = 'results/tables/snora81_overexpression/{models2}_SNORA81_label_{manual_iteration}.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/predict_snora81s_label.py"
