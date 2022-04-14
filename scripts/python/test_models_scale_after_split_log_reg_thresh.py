#!/usr/bin/python3
import pandas as pd
from sklearn import metrics
from sklearn.metrics import roc_curve
from sklearn.linear_model import LogisticRegression
import pickle
import numpy as np

""" Test each model performance on unseen test data and report their accuracy."""
X_train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id_sno')
y_train = pd.read_csv(snakemake.input.y_train, sep='\t')
X_test = pd.read_csv(snakemake.input.X_test[0], sep='\t', index_col='gene_id_sno')
y_test = pd.read_csv(snakemake.input.y_test[0], sep='\t')

# Get best hyperparameters per model
hyperparams_df = pd.read_csv(snakemake.input.best_hyperparameters[0], sep='\t')
hyperparams_df = hyperparams_df.drop('accuracy_cv', axis=1)

def df_to_params(df):
    """ Convert a one-line dataframe into a dict of params and their value. The
        column name corresponds to the key and the value corresponds to
        the value of that param (ex: 'max_depth': 2 where max_depth was the column
        name and 2 was the value in the df)."""
    cols = list(df.columns)
    params = {}
    for col in cols:
        value = df.loc[0, col]
        params[col] = value
    return params

hyperparams = df_to_params(hyperparams_df)




# Define a new class of of LogisticRegression in which we can choose the log_reg threshold used to predict
class LogisticRegressionWithThreshold(LogisticRegression):
    def predict(self, X, threshold=None):
        if threshold == None: # If no threshold passed in, simply call the base class predict, effectively threshold=0.5
            return LogisticRegression.predict(self, X)
        else:
            y_scores = LogisticRegression.predict_proba(self, X)[:, 1]
            y_pred_with_threshold = (y_scores >= threshold).astype(int)

            return y_pred_with_threshold

    def threshold_from_optimal_tpr_minus_fpr(self, X, y):
        # Find optimal log_reg threshold where we maximize the True positive rate (TPR) and minimize the False positive rate (FPR)
        y_scores = LogisticRegression.predict_proba(self, X)[:, 1]
        fpr, tpr, thresholds = roc_curve(y, y_scores)

        optimal_idx = np.argmax(tpr - fpr)

        return thresholds[optimal_idx], tpr[optimal_idx] - fpr[optimal_idx]


# Instantiate that new class and fit the parameters (train on training set)
lrt = LogisticRegressionWithThreshold(C=hyperparams['C'], solver=hyperparams['solver'],
                                random_state=42, max_iter=500, penalty='l2')  # l2 is the default regularization
lrt.fit(X_train, y_train.values.ravel())

# Pickle the model as a .sav file ('wb' for write in binary)
pickle.dump(lrt, open(snakemake.output.pickled_trained_model, 'wb'))

# Find optimal threshold and predict using that threshold instead of 0.5
threshold, optimal_tpr_minus_fpr = lrt.threshold_from_optimal_tpr_minus_fpr(X_test, y_test)
print('Optimal threshold and tpr-fpr:')
print(threshold, optimal_tpr_minus_fpr)
y_pred_thresh = lrt.predict(X_test, threshold)
print('Accuracy:')
print(metrics.accuracy_score(y_test, y_pred_thresh))

# Save training accuracy into df
accu = {}
accu['log_reg_thresh_training_accuracy'] = lrt.score(X_train, y_train.values.ravel())
accu_df = pd.DataFrame(accu, index=[0])
accu_df.to_csv(snakemake.output.training_accuracy, sep='\t', index=False)

# Save test accuracy into df
acc = {}
acc['log_reg_thresh_test_accuracy'] = metrics.accuracy_score(y_test, y_pred_thresh)
acc_df = pd.DataFrame(acc, index=[0])
acc_df.to_csv(snakemake.output.test_accuracy, sep='\t', index=False)
