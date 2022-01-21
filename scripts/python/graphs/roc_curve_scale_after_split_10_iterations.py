#!/usr/bin/python3
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft
import numpy as np

model_colors = snakemake.params.model_colors_dict

# Load test set and labels of all 10 iterations
X_all_test_path, y_all_test_path = snakemake.input.X_test, snakemake.input.y_test
X_test, y_test = [], []
for i, test_path in enumerate(X_all_test_path):
    X_test_iteration = pd.read_csv(test_path, sep='\t', index_col='gene_id_sno')
    y_test_iteration = pd.read_csv(y_all_test_path[i], sep='\t')
    X_test.append(X_test_iteration)
    y_test.append(y_test_iteration)

# Get the name of all models in a list
pickled_models = snakemake.input.pickled_trained_model
model_name = []
for model_path in pickled_models:
    name = model_path.split('_trained')[0]
    name = name.rsplit('/')[-1]
    if name not in model_name:
        model_name.append(name)

# Unpickle and thus instantiate the trained models for all 10 iterations into a dict
loaded_models = {}
for model in model_name:
    iterations_per_model = [path for path in pickled_models if model in path]
    unpickled_models = []
    for iteration_path in iterations_per_model:
        unpickled_model = pickle.load(open(iteration_path, 'rb'))
        unpickled_models.append(unpickled_model)
    loaded_models[model] = unpickled_models


# Compute average true positive rate, false positive rate, AUCs, and stdev of
# AUCs for each model across the 10 iterations
# We use 100 defined thresholds of probabilities to be able to compute a mean and stdev across 10 iterations;
# otherwise, the fpr and tpr given by viz (see below) can have different shape for each iteration
# (because multiple thresholds can have the same x,y --> fpr, tpr coordinates), which would not work to compute an avg per threshold because there would be missing data points
mean_fpr = np.linspace(0, 1, 100)  # we want to create 100 points between 0 and 1 (x-axis) for the roc curve based on
                                    # a smaller number of coordinates given by roc_curve_display (where x,y --> fpr, tpr) using interpolation.
                                    # mean_fpr (mean false positive rate) are 100 evenly distributed x values
mean_tprs = []  # this will contain the avg tpr values (across iterations) per model
std_tprs = []  # this will contain the stdev of tpr values (across iterations) per model
mean_aucs = []  # this will contain the avg auc (across iterations) per model
std_auc = []  # this will contain the stdev of the auc (across iterations) per model
for mod_name, loaded_models in loaded_models.items():
    tprs, aucs = [], []  # true positive rates and AUCs for all iterations of a model
    for i, predictor_per_iteration in enumerate(loaded_models):
        viz = plot_roc_curve(predictor_per_iteration, X_test[i],  # used to access fpr and tps, not to plot the roc curve
                                            y_test[i])
        interpolated_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)  # we interpolate the same number of values (100) for each iteration
        interpolated_tpr[0] = 0.0  # we modify the first value which is 0. to 0.0 on the y-axis
        tprs.append(interpolated_tpr)
        aucs.append(viz.roc_auc)
    mean_tpr_per_model = np.mean(tprs, axis=0)
    mean_tpr_per_model[-1] = 1  # we modify the last value to be exactly 1 on the y axis
    std_tpr_per_model = np.std(tprs, axis=0)
    mean_auc_per_model = auc(mean_fpr, mean_tpr_per_model)
    std_auc_per_model = np.std(aucs)
    mean_tprs.append(mean_tpr_per_model)
    std_tprs.append(std_tpr_per_model)
    mean_aucs.append(mean_auc_per_model)
    std_auc.append(std_auc_per_model)

# Plot the roc curve with error clouds below and above
ft.roc_curve_error_fill(model_name, mean_aucs, mean_fpr, mean_tprs, std_tprs,
                        model_colors, snakemake.output.roc_curve)
