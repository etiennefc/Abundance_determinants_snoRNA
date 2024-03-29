import os


### Split the entire dataset of human snoRNAs in two sets: cv set (10%) and
### training set (90%). The CV set wiill serve to determine the best
### hyperparameters of the models while the training set will serve to fit the
### paramters of the models.

rule fill_na_feature_scale_species_prediction_top4:
    """ Fill NA values across columns in the feature df, fixing them at -5
        because these negative values should not interfere. Then, split dataset
        in 2 sets (cv and train) and after do feature scaling using mean
        normalization to numerical feature columns that have large ranges of
        data value within each of these datasets. Mean normalization corresponds
        to the value substracted by the mean of all values divided by the
        standard deviation of all values((x - mean)/std).The cv and train sets
        are composed of respectively 10 and 90% of all snoRNAs."""
    input:
        feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df
    output:
        cv = 'results/tables/all_features_labels_scaled_cv_top4_species_prediction.tsv',
        train = 'results/tables/all_features_labels_scaled_train_top4_species_prediction.tsv',
        y_cv = 'results/tables/labels_cv_top4_species_prediction.tsv',
        y_train = 'results/tables/labels_train_top4_species_prediction.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/scale_features_species_prediction_top4.py"

rule hyperparameter_tuning_cv_species_prediction_top4:
    """ Tune the hyperparameters of each models (Logistic regression, Support
        vector classifier, Random Forest, KNN and Gradient boosting classifier)
        before even training them, using GridSearchCV with stratified k-fold. The
        CV will take place as a stratifed k-fold (3 fold) using GridSearchCV and
        will return the best hyperparameters for each tested model. """
    input:
        X_cv = rules.fill_na_feature_scale_species_prediction_top4.output.cv,
        y_cv = rules.fill_na_feature_scale_species_prediction_top4.output.y_cv
    output:
        best_hyperparameters = os.path.join(config['path']['hyperparameter_tuning'],
                                    '{models2}_best_params_top4_species_prediction.tsv')
    params:
        hyperparameter_space = lambda wildcards: config['hyperparameter_space'][wildcards.models2]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/hyperparameter_tuning_cv_scale_after_split.py"

rule train_models_species_prediction_top4:
    """ Train (fit) each model on the training set using the best
        hyperparameters found by hyperparameter_tuning_cv_scale_after_split.
        Pickle these fitted models (into .sav files) so that they can be reused
        after without the need to retrain them all over again."""
    input:
        X_train = rules.fill_na_feature_scale_species_prediction_top4.output.train,
        y_train = rules.fill_na_feature_scale_species_prediction_top4.output.y_train,
        best_hyperparameters = rules.hyperparameter_tuning_cv_species_prediction_top4.output.best_hyperparameters
    output:
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    '{models2}_trained_top4_species_prediction.sav'),
        training_accuracy = os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_top4_species_prediction.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/train_models_scale_after_split.py"

rule predict_mouse_snoRNA_label_species_prediction_top4:
    """ Predict the abundance status of all mouse snoRNA based on the
        top4 features (combined_box_hamming. sno_mfe, terminal_stem_mfe and
        host_expressed). We scale the feature_df using the same parameters (mean
        and stdev) used to scale each training set."""
    input:
        feature_df = rules.merge_features_label_mouse.output.feature_df,
        human_snoRNA_feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        pickled_trained_model = rules.train_models_species_prediction_top4.output.pickled_trained_model
    output:
        predicted_label_df = 'results/tables/mouse_prediction/{models2}_predicted_label_species_prediction_top4.tsv',
        scaled_feature_df = 'results/tables/mouse_prediction/{models2}_scaled_features_species_prediction_top4.tsv',
        label_df = 'results/tables/mouse_prediction/{models2}_label_species_prediction_top4.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/predict_species_snoRNA_label.py"

rule test_accuracy_species_prediction_top4:
    """ Test model performance on mouse snoRNA data and return their accuracies."""
    input:
        X_test = rules.predict_mouse_snoRNA_label_species_prediction_top4.output.scaled_feature_df,
        y_test = rules.predict_mouse_snoRNA_label_species_prediction_top4.output.label_df,
        pickled_trained_model = rules.train_models_species_prediction_top4.output.pickled_trained_model
    output:
        test_accuracy = os.path.join(config['path']['test_accuracy_mouse'],
                                    '{models2}_test_accuracy_top4_species_prediction.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/test_models_scale_after_split.py"

rule confusion_matrix_f1_species_prediction_top4:
    """ Create a confusion matrix (all True/False positives or negatives) for
        each model in order to compare whether they misclassify the same snoRNAs
        or not. Compute also the F1 score per model. Do this on mouse snoRNAs."""
    input:
        X_test = rules.predict_mouse_snoRNA_label_species_prediction_top4.output.scaled_feature_df,
        y_test = rules.predict_mouse_snoRNA_label_species_prediction_top4.output.label_df,
        pickled_trained_model = rules.train_models_species_prediction_top4.output.pickled_trained_model
    output:
        confusion_matrix = os.path.join(config['path']['confusion_matrix_f1_mouse'],
                        '{models2}_confusion_matrix_w_f1_score_top4_species_prediction.tsv'),
        info_df = os.path.join(config['path']['confusion_matrix_f1_mouse'],
                        '{models2}_confusion_matrix_values_per_sno_top4_species_prediction.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/confusion_matrix_f1_scale_after_split.py"
