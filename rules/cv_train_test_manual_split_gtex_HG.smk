import os


### Split the entire datasets in 10 stratified folds (10% each). Each of these fold will
### serve as a Test set once (so each snoRNA will be predicted once). The
### other 9 folds will be merged and then split into CV and Training sets (80% and 10%)

### Determine if using host_gene_cutoff from GTEx data instead of from our TGIRT-Seq
### data is the same in terms of accuracy of prediction and feature rank

rule one_hot_encode_before_split_gtex_HG:
    """ One-hot encode categorical features in the feature df (using
        OneHotEncoder) and also label-encode labels BEFORE splitting into
        cv/train/test sets."""
    input:
        feature_df = rules.merge_feature_df_gtex_HG.output.feature_df
    output:
        one_hot_encoded_df = config['path']['one_hot_encoded_df_before_split_gtex_HG']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/one_hot_encode.py"

rule one_hot_encode_before_split_gtex_HG_unpaired:
    """ One-hot encode categorical features in the feature df (using
        OneHotEncoder) and also label-encode labels BEFORE splitting into
        cv/train/test sets."""
    input:
        feature_df = rules.merge_feature_df_gtex_HG_unpaired.output.feature_df
    output:
        one_hot_encoded_df = config['path']['one_hot_encoded_df_before_split_gtex_HG_unpaired']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/one_hot_encode.py"

rule fill_na_feature_scaling_after_manual_split_gtex_HG:
    """ Fill NA values across columns in the feature df, fixing them at -5
        because these negative values should not interfere. Then, split dataset
        in 3 sets (train, cv and test) and after do feature scaling using mean
        normalization to numerical feature columns that have large ranges of
        data value within each of these datasets. Mean normalization corresponds
        to the value substracted by the mean of all values divided by the
        standard deviation of all values((x - mean)/std). Do this to create 10
        subsets of cv-train-test datasets where the combination of all 10 test
        sets includes each snoRNA once. The cv, train and test sets are composed
        of respectively 10, 80 and 10% of all snoRNAs."""
    input:
        feature_df = rules.one_hot_encode_before_split_gtex_HG.output.one_hot_encoded_df
    output:
        cv = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'all_features_labels_scaled_cv_gtex_HG_{manual_iteration}.tsv'),
        train = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'all_features_labels_scaled_train_gtex_HG_{manual_iteration}.tsv'),
        test = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'all_features_labels_scaled_test_gtex_HG_{manual_iteration}.tsv'),
        y_cv = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'labels_cv_gtex_HG_{manual_iteration}.tsv'),
        y_train = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'labels_train_gtex_HG_{manual_iteration}.tsv'),
        y_test = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'labels_test_gtex_HG_{manual_iteration}.tsv')
    params:
        random_state = 42
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/scale_features_after_manual_split.py"

rule hyperparameter_tuning_cv_scale_after_manual_split_gtex_HG:
    """ Tune the hyperparameters of each models (Logistic regression, Support
        vector classifier, Random Forest and Gradient boosting classifier)
        before even training them, using GridSearchCV with stratified k-fold. The
        CV will take place as a stratifed k-fold (3 fold) using GridSearchCV and
        will return the best hyperparameters for each tested model. """
    input:
        X_cv = rules.fill_na_feature_scaling_after_manual_split_gtex_HG.output.cv,
        y_cv = rules.fill_na_feature_scaling_after_manual_split_gtex_HG.output.y_cv
    output:
        best_hyperparameters = os.path.join(config['path']['hyperparameter_tuning'],
                                    '{models2}_best_params_scale_after_split_gtex_HG_{manual_iteration}.tsv')
    params:
        hyperparameter_space = lambda wildcards: config['hyperparameter_space'][wildcards.models2]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/hyperparameter_tuning_cv_scale_after_split.py"

rule train_models_scale_after_manual_split_gtex_HG:
    """ Train (fit) each model on the training set using the best
        hyperparameters found by hyperparameter_tuning_cv_scale_after_split.
        Pickle these fitted models (into .sav files) so that they can be reused
        after without the need to retrain them all over again."""
    input:
        X_train = rules.fill_na_feature_scaling_after_manual_split_gtex_HG.output.train,
        y_train = rules.fill_na_feature_scaling_after_manual_split_gtex_HG.output.y_train,
        best_hyperparameters = rules.hyperparameter_tuning_cv_scale_after_manual_split_gtex_HG.output.best_hyperparameters
    output:
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split_gtex_HG_{manual_iteration}.sav'),
        training_accuracy = os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_scale_after_split_gtex_HG_{manual_iteration}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/train_models_scale_after_split.py"

rule test_models_scale_after_manual_split_gtex_HG:
    """ Test model performance on unseen test data and return their accuracies."""
    input:
        X_test = rules.fill_na_feature_scaling_after_manual_split_gtex_HG.output.test,
        y_test = rules.fill_na_feature_scaling_after_manual_split_gtex_HG.output.y_test,
        pickled_trained_model = rules.train_models_scale_after_manual_split_gtex_HG.output.pickled_trained_model
    output:
        test_accuracy = os.path.join(config['path']['test_accuracy'],
                                    '{models2}_test_accuracy_scale_after_split_gtex_HG_{manual_iteration}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/test_models_scale_after_split.py"

rule confusion_matrix_f1_scale_after_manual_split_gtex_HG:
    """ Create a confusion matrix (all True/False positives or negatives) for
        each model in order to compare whether they misclassify the same snoRNAs
        or not. Compute also the F1 score per model."""
    input:
        X_test = rules.fill_na_feature_scaling_after_manual_split_gtex_HG.output.test,
        y_test = rules.fill_na_feature_scaling_after_manual_split_gtex_HG.output.y_test,
        pickled_trained_model = rules.train_models_scale_after_manual_split_gtex_HG.output.pickled_trained_model
    output:
        confusion_matrix = os.path.join(config['path']['confusion_matrix_f1'],
                        '{models2}_confusion_matrix_w_f1_score_scale_after_split_gtex_HG_{manual_iteration}.tsv'),
        info_df = os.path.join(config['path']['confusion_matrix_f1'],
                        '{models2}_confusion_matrix_values_per_sno_scale_after_split_gtex_HG_{manual_iteration}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/confusion_matrix_f1_scale_after_split.py"
