#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)  # ignore all user warnings
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import functions as ft

X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')
y_test = pd.read_csv(snakemake.input.y_test, sep='\t')

# Unpickle and thus instantiate the 4 trained models
log_reg = pickle.load(open(snakemake.input.log_reg, 'rb'))
svc = pickle.load(open(snakemake.input.svc, 'rb'))
rf = pickle.load(open(snakemake.input.rf, 'rb'))
gbm = pickle.load(open(snakemake.input.gbm, 'rb'))
knn = pickle.load(open(snakemake.input.knn, 'rb'))
# Create the ROC curve
classifiers = [log_reg, svc, rf, gbm, knn]
ft.roc_curve(classifiers, X_test, y_test, "False positive rate",
            "True positive rate", "ROC curves of the 5 models showing"+"\n"+"their performance on the test dataset",
            snakemake.output.roc_curve)
