#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft
import re
import statistics as st

# Define function to return the average and standard deviation of accuracies of the 20 iterations per model in a respective dict each
def get_avg_stdev(dict_of_all_models_accuracies):
    avg, stdev = {}, {}
    for i in range(0, len(sorted(dict_of_all_models_accuracies.keys())), 20):  # sort to regroup in order all 20 iterations per model
        iterations_per_model = sorted(dict_of_all_models_accuracies.keys())[i:i+20]  # select the 20 iterations names per model
        accuracies_per_model = [dict_of_all_models_accuracies[iteration] for iteration in iterations_per_model]  # select corresponding accuracies of these 20 iterations
        model_name = iterations_per_model[0].split('_')[0]  # Get the model name
        avg_acc, stdev_acc = st.mean(accuracies_per_model), st.stdev(accuracies_per_model)
        avg[model_name], stdev[model_name] = avg_acc, stdev_acc
    return avg, stdev

# Get the accuracy of all models on the CV set
cv_accuracy = {}
for i, model in enumerate(snakemake.input.cv_accuracy):
    df = pd.read_csv(model, sep='\t')
    substring_model = re.findall("(log_reg|svc|rf|gbm|knn)", model)[0]
    iteration = model.split('_')[-2]
    accuracy = float(df['accuracy_cv'].values)
    cv_accuracy[f'{substring_model}_{iteration}'] = accuracy

# Get the accuracy of all models on the training set
train_accuracy = {}
for i, model in enumerate(snakemake.input.training_accuracy):
    df = pd.read_csv(model, sep='\t')
    substring_model = re.findall("(log_reg|svc|rf|gbm|knn)", model)[0]
    iteration = model.split('_')[-2]
    accuracy = float(df[f'{substring_model}_training_accuracy'].values)
    train_accuracy[f'{substring_model}_{iteration}'] = accuracy

# Get the accuracy of all models on the test set
test_accuracy = {}
for i, model in enumerate(snakemake.input.test_accuracy):
    df = pd.read_csv(model, sep='\t')
    substring_model = re.findall("(log_reg|svc|rf|gbm|knn)", model)[0]
    iteration = model.split('_')[-2]
    accuracy = float(df[f'{substring_model}_test_accuracy'].values)
    test_accuracy[f'{substring_model}_{iteration}'] = accuracy

# Get the average and standard deviation of CV, training and test sets accuracies for all models across the 20 iterations
cv_avg, cv_stdev = get_avg_stdev(cv_accuracy)
train_avg, train_stdev = get_avg_stdev(train_accuracy)
test_avg, test_stdev = get_avg_stdev(test_accuracy)

print(cv_accuracy)
print(train_accuracy)
print(test_accuracy)
# Create the df of all accuracies
all_accuracies = pd.DataFrame.from_dict([cv_avg, train_avg, test_avg])
all_accuracies.index = ["cv", "train", "test"]
all_accuracies = all_accuracies.transpose()
all_accuracies_hue = all_accuracies.copy()
all_accuracies_hue['model'] = all_accuracies_hue.index

# Create the connected scatter plot
color_dict = snakemake.params.colors
color_dict['log'] = color_dict.pop('log_reg')  # patch to link the good color to log_reg
ft.connected_scatter_errbars(all_accuracies, all_accuracies_hue, 'model',
                    color_dict, 'Dataset', [cv_stdev, train_stdev, test_stdev],
                    ['Cross-\nvalidation', 'Training', 'Test'], 'Dataset', 'Accuracy',
                    snakemake.output.scatter)
