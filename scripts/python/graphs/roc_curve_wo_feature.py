#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
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
X_train, X_test, y_train, y_test = train_test_split(X_total_train, y_total_train,
                                    test_size=232, train_size=1077, random_state=42,
                                    stratify=y_total_train)


# Unpickle and thus instantiate the 4 trained models
log_reg = pickle.load(open(snakemake.input.log_reg, 'rb'))
svc = pickle.load(open(snakemake.input.svc, 'rb'))
rf = pickle.load(open(snakemake.input.rf, 'rb'))
gbm = pickle.load(open(snakemake.input.gbm, 'rb'))
knn = pickle.load(open(snakemake.input.knn, 'rb'))
# Create the ROC curve
classifiers = [log_reg, svc, rf, gbm, knn]
ft.roc_curve(classifiers, X_test, y_test, "False positive rate",
            "True positive rate", "ROC curves of the 5 models showing their performance on"+"\n"+"the test dataset without the "+snakemake.wildcards.feature_effect+" feature",
            snakemake.output.roc_curve)
