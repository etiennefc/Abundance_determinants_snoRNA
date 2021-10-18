import os

include: "feature_normalization.smk"
###CV, train and test models based on the datasets generated with the scaling after the split of the 3 datasets

rule hyperparameter_tuning_cv_scale_after_split:
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
        X_cv = rules.fill_na_feature_scaling_after_split.output.cv,
        y_cv = rules.fill_na_feature_scaling_after_split.output.y_cv
    output:
        best_hyperparameters = os.path.join(config['path']['hyperparameter_tuning'],
                                    '{models2}_best_params_scale_after_split.tsv')
    params:
        hyperparameter_space = lambda wildcards: config['hyperparameter_space'][wildcards.models2]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/hyperparameter_tuning_cv_scale_after_split.py"

rule train_models_scale_after_split:
    """ Train (fit) each model on the training set using the best
        hyperparameters found by hyperparameter_tuning_cv_scale_after_split.
        Pickle these fitted models (into .sav files) so that they can be reused
        after without the need to retrain them all over again."""
    input:
        X_train = rules.fill_na_feature_scaling_after_split.output.train,
        y_train = rules.fill_na_feature_scaling_after_split.output.y_train,
        best_hyperparameters = rules.hyperparameter_tuning_cv_scale_after_split.output.best_hyperparameters
    output:
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav'),
        training_accuracy = os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_scale_after_split.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/train_models_scale_after_split.py"

rule test_models_scale_after_split:
    """ Test model performance on unseen test data and return their accuracies."""
    input:
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        y_test = rules.fill_na_feature_scaling_after_split.output.y_test,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    output:
        test_accuracy = os.path.join(config['path']['test_accuracy'],
                                    '{models2}_test_accuracy_scale_after_split.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/test_models_scale_after_split.py"

rule confusion_matrix_f1_scale_after_split:
    """ Create a confusion matrix (all True/False positives or negatives) for
        each model in order to compare whether they misclassify the same snoRNAs
        or not. Compute also the F1 score per model."""
    input:
        X_test = rules.fill_na_feature_scaling_after_split.output.test,
        y_test = rules.fill_na_feature_scaling_after_split.output.y_test,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_scale_after_split.sav')
    output:
        confusion_matrix = os.path.join(config['path']['confusion_matrix_f1'],
                        '{models2}_confusion_matrix_w_f1_score_scale_after_split.tsv'),
        info_df = os.path.join(config['path']['confusion_matrix_f1'],
                        '{models2}_confusion_matrix_values_per_sno_scale_after_split.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/confusion_matrix_f1_scale_after_split.py"
