#!/usr/bin/python3
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import metrics
import pickle

""" Test each model performance on unseen test data and report their accuracy."""

X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')
y_test = pd.read_csv(snakemake.input.y_test, sep='\t')

# Unpickle and thus instantiate the trained model defined by the 'models' wildcard
model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))

# Predict label (expressed (1) or not_expressed (1)) on test data and compare to y_test
y_pred = model.predict(X_test)
print(snakemake.wildcards.models2)
print(metrics.accuracy_score(y_test, y_pred))
acc = {}
acc[snakemake.wildcards.models2+'_test_accuracy'] = metrics.accuracy_score(y_test, y_pred)
acc_df = pd.DataFrame(acc, index=[0])
acc_df.to_csv(snakemake.output.test_accuracy, sep='\t', index=False)
