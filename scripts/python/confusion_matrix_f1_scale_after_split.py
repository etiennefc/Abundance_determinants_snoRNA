#!/usr/bin/python3
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, f1_score
import pickle

""" Return the confusion matrix associated to each model and their F1 score.
    Return also a df containing each snoRNA in the test set and their associated
    confusion matrix value (i.e. TN, FP, FN or TP)."""

X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id_sno')
y_test = pd.read_csv(snakemake.input.y_test, sep='\t')
y_test.index = X_test.index  # set gene_id_sno as index

# Unpickle and thus instantiate the trained model defined by the 'models' wildcard
model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))

# Predict label (expressed (1) or not_expressed (0)) on test data and compare to y_test
y_pred = model.predict(X_test)

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

info_df.loc[(info_df['label'] == 0) & (info_df['predicted_label'] == 0), 'confusion_matrix_val_' + snakemake.wildcards.models2] = 'TN'
info_df.loc[(info_df['label'] == 0) & (info_df['predicted_label'] == 1), 'confusion_matrix_val_' + snakemake.wildcards.models2] = 'FP'
info_df.loc[(info_df['label'] == 1) & (info_df['predicted_label'] == 0), 'confusion_matrix_val_' + snakemake.wildcards.models2] = 'FN'
info_df.loc[(info_df['label'] == 1) & (info_df['predicted_label'] == 1), 'confusion_matrix_val_' + snakemake.wildcards.models2] = 'TP'

info_df.to_csv(snakemake.output.info_df, sep='\t')
