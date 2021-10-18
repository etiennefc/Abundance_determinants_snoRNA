#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import shap
import numpy as np
# Generate the same CV, training and test sets (only the test set will be
# used in this script) that were generated in hyperparameter_tuning_cv and train_models
# (respectively 15%, 70% and 15% of all dataset examples)
df = pd.read_csv(snakemake.input.df, sep='\t', index_col='gene_id_sno')
X = df.drop('label', axis=1)
y = df['label']

# First the CV vs total_train split
X_total_train, X_cv, y_total_train, y_cv = train_test_split(X, y, test_size=0.15,
                                            random_state=42, stratify=y)

# Next the total_train is split into train and test sets (1077 and 232 correspond
# to the number of examples in train and test sets respectively to get an
# approximately 70 % and 15 % of all examples in these two datasets)
X_train, X_test, y_train, y_test = train_test_split(X_total_train, y_total_train, test_size=232, train_size=1077, random_state=42, stratify=y_total_train)


# Unpickle and thus instantiate the model represented by the 'models' wildcard
# Instantiate the explainer using the X_train as background data and X_test to generate shap local values for one snoRNA
if snakemake.wildcards.models == "log_reg":
    model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    explainer = shap.LinearExplainer(model, shap.sample(X_train, 100, random_state=42))  # reduce number of background sample to 100
    shap_values = explainer.shap_values(X_test)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    #shap.plots.bar(shap_values, show=False, max_display=50)
    shap.summary_plot(shap_values, X_test, plot_type='bar', max_display=50, show=False)
    plt.savefig(snakemake.output.bar_plot, bbox_inches='tight', dpi=600)

else:
    model2 = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    explainer2 = shap.KernelExplainer(model2.predict, shap.sample(X_train, 100, random_state=42)) # reduce number of background sample to 100
    shap_values2 = explainer2.shap_values(X_test)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    #shap.plots.bar(shap_values2, show=False, max_display=50)
    shap.summary_plot(shap_values2, X_test, plot_type='bar', max_display=50, show=False)
    plt.savefig(snakemake.output.bar_plot, bbox_inches='tight', dpi=600)
