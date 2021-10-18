import os

include: "feature_normalization.smk"

### Redo CV, train and test, but without the one or all of the top 3 most predictive features
### to see how it affects the performance of the models

rule remove_feature:
    """ Remove one or multiple features from the original dataset."""
    input:
        df = rules.one_hot_encode.output.one_hot_encoded_df
    output:
        df_wo_feature = os.path.join(
                            config['path']['one_hot_encoded_df_wo_feature'],
                            '{feature_effect}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/remove_feature.py"

rule hyperparameter_tuning_cv_wo_feature:
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
        df = rules.remove_feature.output.df_wo_feature
    output:
        best_hyperparameters = os.path.join(config['path']['hyperparameter_tuning'],
                                    '{models2}_best_params_wo_{feature_effect}.tsv')
    params:
        hyperparameter_space = lambda wildcards: config['hyperparameter_space'][wildcards.models2]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/hyperparameter_tuning_cv.py"

rule train_models_wo_feature:
    """ Train (fit) each model on the training set using the best
        hyperparameters found by hyperparameter_tuning_cv. Pickle these fitted
        models (into .sav files) so that they can be reused after without the
        need to retrain them all over again."""
    input:
        df = rules.remove_feature.output.df_wo_feature,
        best_hyperparameters = rules.hyperparameter_tuning_cv_wo_feature.output.best_hyperparameters
    output:
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_wo_{feature_effect}.sav'),
        training_accuracy = os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_wo_{feature_effect}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/train_models.py"

rule test_models_wo_feature:
    """ Test model performance on unseen test data and return their accuracies."""
    input:
        df = rules.remove_feature.output.df_wo_feature,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_wo_{feature_effect}.sav')
    output:
        test_accuracy = os.path.join(config['path']['test_accuracy'],
                                    '{models2}_test_accuracy_wo_{feature_effect}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/test_models.py"

rule confusion_matrix_f1_wo_feature:
    """ Create a confusion matrix (all True/False positives or negatives) for
        each model in order to compare whether they misclassify the same snoRNAs
        or not. Compute also the F1 score per model."""
    input:
        df = rules.remove_feature.output.df_wo_feature,
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models}_trained_wo_{feature_effect}.sav')
    output:
        confusion_matrix = os.path.join(config['path']['confusion_matrix_f1'],
                        '{models}_confusion_matrix_w_f1_score_wo_{feature_effect}.tsv'),
        info_df = os.path.join(config['path']['confusion_matrix_f1'],
                        '{models}_confusion_matrix_values_per_sno_wo_{feature_effect}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/confusion_matrix_f1.py"
