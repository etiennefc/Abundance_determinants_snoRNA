import os



rule fill_na_feature_scaling:
    """ Fill NA values across columns in the feature df, fixing them at -5
        because these negative values should not interfere. Then do feature
        scaling using mean normalization to numerical feature columns that have large
        ranges of data value. Mean normalization corresponds to the value
        substracted by the mean of all values divided by the standard deviation
        of all values((x - mean)/std)."""
    input:
        feature_df = config['path']['feature_df']
    output:
        scaled_feature_df = config['path']['scaled_feature_df']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/scale_features.py"


rule one_hot_encode:
    """ One-hot encode categorical features in the feature df (using
        OneHotEncoder) and also label-encode labels."""
    input:
        feature_df = rules.fill_na_feature_scaling.output.scaled_feature_df
    output:
        one_hot_encoded_df = config['path']['one_hot_encoded_df']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/one_hot_encode.py"


rule one_hot_encode_before_split:
    """ One-hot encode categorical features in the feature df (using
        OneHotEncoder) and also label-encode labels."""
    input:
        feature_df = config['path']['feature_df']
    output:
        one_hot_encoded_df = config['path']['one_hot_encoded_df_before_split']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/one_hot_encode.py"

rule fill_na_feature_scaling_after_split:
    """ Fill NA values across columns in the feature df, fixing them at -5
        because these negative values should not interfere. Then, split dataset
        in 3 sets (train, cv and test) and after do feature scaling using mean
        normalization to numerical feature columns that have large ranges of
        data value within each of these datasets. Mean normalization corresponds
        to the value substracted by the mean of all values divided by the
        standard deviation of all values((x - mean)/std)."""
    input:
        feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df
    output:
        cv = config['path']['scaled_feature_cv_df'],
        train = config['path']['scaled_feature_train_df'],
        test = config['path']['scaled_feature_test_df'],
        y_cv = config['path']['scaled_feature_y_cv'],
        y_train = config['path']['scaled_feature_y_train'],
        y_test = config['path']['scaled_feature_y_test']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/scale_features_after_split.py"
