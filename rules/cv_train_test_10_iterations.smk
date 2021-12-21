import os


###CV, train and test models based on the datasets generated with the scaling after the split of the 3 datasets
# Do this cv-train-test 10 times (iterations) to see if models are comparable depending on the initial split of the 3 datasets

rule fill_na_feature_scaling_after_split_10_iterations:
    """ Fill NA values across columns in the feature df, fixing them at -5
        because these negative values should not interfere. Then, split dataset
        in 3 sets (train, cv and test) and after do feature scaling using mean
        normalization to numerical feature columns that have large ranges of
        data value within each of these datasets. Mean normalization corresponds
        to the value substracted by the mean of all values divided by the
        standard deviation of all values((x - mean)/std). Do this to create 10
        subsets of cv-train-test datasets."""
    input:
        feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df
    output:
        cv = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'all_features_labels_scaled_cv_{iteration}.tsv'),
        train = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'all_features_labels_scaled_train_{iteration}.tsv'),
        test = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'all_features_labels_scaled_test_{iteration}.tsv'),
        y_cv = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'labels_cv_{iteration}.tsv'),
        y_train = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'labels_train_{iteration}.tsv'),
        y_test = os.path.join(config['path']['scaled_feature_10_iterations'],
                        'labels_test_{iteration}.tsv')
    params:
        random_state = config['random_state_iteration']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/scale_features_after_split_10_iterations.py"

rule hyperparameter_tuning_cv_scale_after_split_10_iterations:
    """ Tune the hyperparameters of each models (Logistic regression, Support
        vector classifier, Random Forest and Gradient boosting classifier)
        before even training them, using GridSearchCV with stratified k-fold.
        First, shuffle the dataset and split it into Cross-validation (CV) set
        and training set (respectively 15 % and 85 % of all examples). The CV
        will take place as a stratifed k-fold (3 fold) using GridSearchCV and
        will return the best hyperparameters for each tested model. Of note, the
        training set will be also split later on into training and test set
        (respectively 60 % and 15 % of all examples)."""
    input:
        X_cv = rules.fill_na_feature_scaling_after_split_10_iterations.output.cv,
        y_cv = rules.fill_na_feature_scaling_after_split_10_iterations.output.y_cv
    output:
        best_hyperparameters = os.path.join(config['path']['hyperparameter_tuning'],
                                    '{models2}_best_params_scale_after_split_{iteration}.tsv')
    params:
        hyperparameter_space = lambda wildcards: config['hyperparameter_space'][wildcards.models2]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/hyperparameter_tuning_cv_scale_after_split.py"

rule train_models_scale_after_split_10_iterations:
    """ Train (fit) each model on the training set using the best
        hyperparameters found by hyperparameter_tuning_cv_scale_after_split.
        Pickle these fitted models (into .sav files) so that they can be reused
        after without the need to retrain them all over again."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split_10_iterations.output.train,
        y_train = rules.fill_na_feature_scaling_after_split_10_iterations.output.y_train,
        best_hyperparameters = rules.hyperparameter_tuning_cv_scale_after_split_10_iterations.output.best_hyperparameters
    output:
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split_{iteration}.sav'),
        training_accuracy = os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_scale_after_split_{iteration}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/train_models_scale_after_split.py"

rule test_models_scale_after_split_10_iterations:
    """ Test model performance on unseen test data and return their accuracies."""
    input:
        X_test = rules.fill_na_feature_scaling_after_split_10_iterations.output.test,
        y_test = rules.fill_na_feature_scaling_after_split_10_iterations.output.y_test,
        pickled_trained_model = rules.train_models_scale_after_split_10_iterations.output.pickled_trained_model
    output:
        test_accuracy = os.path.join(config['path']['test_accuracy'],
                                    '{models2}_test_accuracy_scale_after_split_{iteration}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/test_models_scale_after_split.py"

rule confusion_matrix_f1_scale_after_split_10_iterations:
    """ Create a confusion matrix (all True/False positives or negatives) for
        each model in order to compare whether they misclassify the same snoRNAs
        or not. Compute also the F1 score per model."""
    input:
        X_test = rules.fill_na_feature_scaling_after_split_10_iterations.output.test,
        y_test = rules.fill_na_feature_scaling_after_split_10_iterations.output.y_test,
        pickled_trained_model = rules.train_models_scale_after_split_10_iterations.output.pickled_trained_model
    output:
        confusion_matrix = os.path.join(config['path']['confusion_matrix_f1'],
                        '{models2}_confusion_matrix_w_f1_score_scale_after_split_{iteration}.tsv'),
        info_df = os.path.join(config['path']['confusion_matrix_f1'],
                        '{models2}_confusion_matrix_values_per_sno_scale_after_split_{iteration}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/confusion_matrix_f1_scale_after_split.py"
