import os


### Split the entire dataset of human snoRNAs in two sets: cv set (10%) and
### training set (90%). The CV set will serve to determine the best
### hyperparameters of the log_reg model while the training set will serve to fit the
### parameters of the models.

rule train_test_accuracy_species_prediction_top4_log_reg_thresh:
    """ Train log_reg model and test using the optimal log_reg threshold to
        maximize true positive rate while minimizing the false positive rate.
        Then test model performance on mouse snoRNA data and return its accuracy."""
    input:
        best_hyperparameters = expand(rules.hyperparameter_tuning_cv_species_prediction_top4_rs.output.best_hyperparameters, models2=['log_reg'], allow_missing=True),
        X_train = rules.fill_na_feature_scale_species_prediction_top4_rs.output.train,
        y_train = rules.fill_na_feature_scale_species_prediction_top4_rs.output.y_train,
        X_test = expand(rules.predict_mouse_snoRNA_label_species_prediction_top4.output.scaled_feature_df, models2=['log_reg']),
        y_test = expand(rules.predict_mouse_snoRNA_label_species_prediction_top4.output.label_df, models2=['log_reg'])
    output:
        training_accuracy = os.path.join(config['path']['training_accuracy'],
                                    'log_reg_training_accuracy_top4_species_prediction_thresh_{rs}.tsv'),
        test_accuracy = os.path.join(config['path']['test_accuracy_mouse'],
                                    'log_reg_test_accuracy_top4_species_prediction_thresh_{rs}.tsv'),
        pickled_trained_model = os.path.join(config['path']['trained_models'],
                                    'log_reg_trained_top4_species_prediction_thresh_{rs}.sav')

    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/test_models_scale_after_split_log_reg_thresh.py"

rule confusion_matrix_f1_species_prediction_top4_log_reg_thresh:
    """ Create a confusion matrix (all True/False positives or negatives) for
        each model in order to compare whether they misclassify the same snoRNAs
        or not. Compute also the F1 score per model. Do this on mouse snoRNAs
        for the log_reg thresh model."""
    input:
        X_test = expand(rules.predict_mouse_snoRNA_label_species_prediction_top4.output.scaled_feature_df, models2=['log_reg']),
        y_test = expand(rules.predict_mouse_snoRNA_label_species_prediction_top4.output.label_df, models2=['log_reg']),
        pickled_trained_model = rules.train_test_accuracy_species_prediction_top4_log_reg_thresh.output.pickled_trained_model
    output:
        confusion_matrix = os.path.join(config['path']['confusion_matrix_f1_mouse'],
                        'log_reg_confusion_matrix_w_f1_score_top4_species_prediction_thresh_{rs}.tsv'),
        info_df = os.path.join(config['path']['confusion_matrix_f1_mouse'],
                        'log_reg_confusion_matrix_values_per_sno_top4_species_prediction_thresh_{rs}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/confusion_matrix_f1_scale_after_split_log_reg_thresh.py"
