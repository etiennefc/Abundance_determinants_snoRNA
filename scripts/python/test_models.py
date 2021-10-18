#!/usr/bin/python3
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import metrics
import pickle

""" Test each model performance on unseen test data and report their accuracy."""

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
