#!/usr/bin/python3
import pandas as pd
from sklearn.metrics import confusion_matrix, f1_score, roc_curve, accuracy_score
from sklearn.linear_model import LogisticRegression
import pickle
import numpy as np

""" Return the confusion matrix associated to each model and their F1 score.
    Return also a df containing each snoRNA in the test set and their associated
    confusion matrix value (i.e. TN, FP, FN or TP)."""

X_test = pd.read_csv(snakemake.input.X_test[0], sep='\t', index_col='gene_id_sno')
y_test = pd.read_csv(snakemake.input.y_test[0], sep='\t')
y_test.index = X_test.index  # set gene_id_sno as index

# Define class LogisticRegressionWithThreshold
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


# Unpickle and thus instantiate the trained model defined by the 'models' wildcard
model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))

# Find optimal threshold and predict using that threshold instead of 0.5
threshold, optimal_tpr_minus_fpr = model.threshold_from_optimal_tpr_minus_fpr(X_test, y_test)
print('Optimal threshold and tpr-fpr:')
print(threshold, optimal_tpr_minus_fpr)
y_pred = model.predict(X_test, threshold)
print('Accuracy:')
print(accuracy_score(y_test, y_pred))



# Compute the confusion matrix (where T:True, F:False, P:Positive, N:Negative)
TN, FP, FN, TP = confusion_matrix(y_test, y_pred).ravel()

# Compute F1 score
f1 = f1_score(y_test, y_pred)

# Return confusion matrix with F1 socre included
matrix_dict = {'true_negatives': TN, 'false_positives': FP,
                'false_negatives': FN, 'true_positives': TP, 'f1_score': f1}
matrix_df = pd.DataFrame(matrix_dict, index=[0])
matrix_df.to_csv(snakemake.output.confusion_matrix, sep='\t', index=False)

# Return snoRNAs and their confusion matrix value (TN, FP, FN and TP) as a df
y_test_df = pd.DataFrame(y_test)
y_test_df = y_test_df.reset_index()
y_pred_df = pd.DataFrame(y_pred)
y_pred_df.columns = ['predicted_label']

info_df = pd.concat([y_test_df, y_pred_df], axis=1)

info_df.loc[(info_df['label'] == 0) & (info_df['predicted_label'] == 0), 'confusion_matrix_val_log_reg_thresh'] = 'TN'
info_df.loc[(info_df['label'] == 0) & (info_df['predicted_label'] == 1), 'confusion_matrix_val_log_reg_thresh'] = 'FP'
info_df.loc[(info_df['label'] == 1) & (info_df['predicted_label'] == 0), 'confusion_matrix_val_log_reg_thresh'] = 'FN'
info_df.loc[(info_df['label'] == 1) & (info_df['predicted_label'] == 1), 'confusion_matrix_val_log_reg_thresh'] = 'TP'

info_df.to_csv(snakemake.output.info_df, sep='\t')
