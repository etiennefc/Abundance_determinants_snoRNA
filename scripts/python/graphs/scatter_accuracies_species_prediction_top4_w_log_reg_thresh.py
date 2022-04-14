#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft
import re

train_accuracy_df_paths = list(snakemake.input.training_accuracy) + [snakemake.input.training_accuracy_log_reg_thresh[0]]
test_accuracy_df_paths = list(snakemake.input.test_accuracy) + [snakemake.input.test_accuracy_log_reg_thresh[0]]

# Get the accuracy of all models on the CV set
cv_accuracy = {}
for i, model in enumerate(snakemake.input.cv_accuracy):
    df = pd.read_csv(model, sep='\t')
    substring_model = re.findall("(log_reg|svc|rf|gbm|knn)", model)[0]
    accuracy = float(df['accuracy_cv'].values)
    cv_accuracy[substring_model] = accuracy

# Get the accuracy of all models on the training set
train_accuracy = {}
for i, model in enumerate(train_accuracy_df_paths):
    df = pd.read_csv(model, sep='\t')
    substring_model = re.findall("(log_reg|svc|rf|gbm|knn)", model)[0]
    if substring_model == 'log_reg':
        accuracy = float(df[f'{substring_model}_thresh_training_accuracy'].values)
    else:
        accuracy = float(df[f'{substring_model}_training_accuracy'].values)
    train_accuracy[substring_model] = accuracy

# Get the accuracy of all models on the test set
test_accuracy = {}
for i, model in enumerate(test_accuracy_df_paths):
    df = pd.read_csv(model, sep='\t')
    substring_model = re.findall("(log_reg|svc|rf|gbm|knn)", model)[0]
    if substring_model == 'log_reg':
        accuracy = float(df[f'{substring_model}_thresh_test_accuracy'].values)
    else:
        accuracy = float(df[f'{substring_model}_test_accuracy'].values)
    test_accuracy[substring_model] = accuracy


# Create the df of all accuracies
all_accuracies = pd.DataFrame.from_dict([cv_accuracy, train_accuracy, test_accuracy])
all_accuracies.index = ["cv", "train", "test"]
all_accuracies = all_accuracies.transpose()
all_accuracies_hue = all_accuracies.copy()
all_accuracies_hue['model'] = all_accuracies_hue.index

# Create the connected scatter plot
ft.connected_scatter(all_accuracies, all_accuracies_hue, 'model',
                    snakemake.params.colors.values(), 'Dataset',
                    ['Cross-validation', 'Training', 'Test'], 'Dataset', 'Accuracy',
                    snakemake.output.scatter)
