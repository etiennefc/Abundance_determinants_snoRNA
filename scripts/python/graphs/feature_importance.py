#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import numpy as np
import functions as ft

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


# Function to order feature by importance absolute values
def abs_feature_imp(feature_name_list, feature_val_list):
    d = dict(zip(feature_name_list, feature_val_list))
    ordered_val = sorted(feature_val_list, reverse=True, key=abs)  # absolute value used for sorting
    keys = list(d.keys())
    vals = list(d.values())

    ordered_keys = []
    temp_index = 0
    for val in ordered_val:
        if len(list(np.where(ordered_val == val)[0])) > 1:  # if multiple features have the same value
            multiple_keys = [k for k,v in d.items() if v == val]
            key = multiple_keys[temp_index]
            temp_index += 1  # update to select the next feature name with same value
            ordered_keys.append(key)
        else:
            key = keys[vals.index(val)]
            ordered_keys.append(key)

    return ordered_keys, ordered_val


# Unpickle and thus instantiate the model represented by the 'models' wildcard
if snakemake.wildcards.models == "log_reg":
    model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    importance = model.coef_[0]
    x_tick_labels, values = abs_feature_imp(list(X_train.columns), importance)
    ft.simple_bar(values, x_tick_labels, '', 'Feature name',
        'Feature importance', snakemake.output.bar_plot)


elif snakemake.wildcards.models == "svc":
    # Create empty figure, since we can't extract feature importance with an SVC with a sigmoid kernel
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    plt.savefig(snakemake.output.bar_plot, bbox_inches='tight', dpi=600)

else:  # for gbm and rf
    model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))
    importance = model.feature_importances_
    x_tick_labels, values = abs_feature_imp(list(X_train.columns), importance)
    ft.simple_bar(values, x_tick_labels, '', 'Feature name',
        'Feature importance', snakemake.output.bar_plot)
