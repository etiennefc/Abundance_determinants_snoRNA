#!/usr/bin/python3
import pandas as pd
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.utils import shuffle
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder

""" Scale features based on the training sets used to train the models on human
    snoRNAs and predict the abundance status of species snoRNAs."""

random_state = snakemake.params.random_state
human_feature_df_0 = pd.read_csv(snakemake.input.human_snoRNA_feature_df, sep='\t')
feature_df = pd.read_csv(snakemake.input.feature_df, sep='\t')  # species snoRNAs
feature_df_copy = feature_df.copy()
threshold = pd.read_csv(snakemake.input.threshold[0], sep='\t')  # threshold used for log_reg decision
thresh = threshold.iloc[0, 0]

## Species dataset processing
# Drop gene_name col and select numerical feature columns
feature_df = feature_df.drop(columns='gene_name')
numerical_features_label = feature_df.copy()
numerical_features_label = numerical_features_label.select_dtypes(include=['int64', 'float64'])

# Convert categorical features (astype 'object') into numpy array, convert the
# strings in that array into numbers with LabelEncoder and then one-hot encode this numerical array
dfs = [feature_df[['gene_id_sno']]]
feature_df = feature_df.drop(columns=['gene_id_sno'])
df_cat = feature_df.select_dtypes(include=['object'])
cols = df_cat.columns
for i, col in enumerate(cols):
    df_cat = feature_df[[col]]

    # Convert column into numpy array
    array_cat = df_cat.values.reshape(-1, 1)  # -1 infers the length of df_cat

    # Transform string array into numerical array
    le = LabelEncoder()
    df_cat[col+'_cat'] = le.fit_transform(df_cat[col])

    # Get the string that is linked to each numerical category created by LabelEncoder
    label_dict = dict(zip(le.classes_, le.transform(le.classes_)))
    labels = list(label_dict.keys())

    # One-hot encode the numerical array that was created
    enc = OneHotEncoder(handle_unknown='ignore')
    one_hot_array = enc.fit_transform(df_cat[[col+'_cat']]).toarray()
    enc_df = pd.DataFrame(one_hot_array, columns=labels)

    dfs.append(enc_df)

# Concat all one-hot encoded categorical columns
final_df = pd.concat(dfs, axis=1)

# Concat numerical features and label at the end and set index as sno_id
final_df = pd.concat([final_df, numerical_features_label], axis=1)
final_df = final_df.set_index('gene_id_sno')

# Keep only relevant columns
final_df = final_df[['sno_mfe', 'terminal_stem_mfe', 'combined_box_hamming', 'host_expressed']]





## Human dataset processing
# Fill NaN with -5
human_feature_df_0 = human_feature_df_0.fillna(-5)
human_feature_df = human_feature_df_0.copy()
human_feature_df = shuffle(human_feature_df, random_state=random_state)

# Keep only top 4 features
X = human_feature_df.drop('label', axis=1)
X = X[['gene_id_sno', 'sno_mfe', 'terminal_stem_mfe', 'combined_box_hamming', 'host_expressed']]  # same order as in final_df that we want to predict
y = human_feature_df['label']


# Split the X dataset into cv and train test (respectively 10 and 90 % of the total df)
X_train, X_cv, y_train, y_cv = train_test_split(X, y,
                            test_size=0.1, train_size=0.9, random_state=random_state,
                            stratify=y)

# Create the scaler that will be used to mean normalize species snoRNA features
df = X_train.set_index('gene_id_sno')
scaler = StandardScaler().fit(df)



# Scale final_df (all species snoRNAs with top4 features) by the same scaling factor used to scale the training set in human snoRNAs
scaled_feature_df = scaler.transform(final_df)


# Define a new class of LogisticRegression in which we can choose the log_reg threshold used to predict
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


# Unpickle and thus instantiate the trained log_reg thresh model
model = pickle.load(open(snakemake.input.pickled_trained_model[0], 'rb'))

# Predict label (expressed (1) or not_expressed (0)) on species snoRNAs
y_pred = model.predict(scaled_feature_df, thresh)
prediction_df = pd.DataFrame(y_pred, index=final_df.index, columns=['predicted_label'])
prediction_df = prediction_df.reset_index()
prediction_df = prediction_df.replace(to_replace=[0, 1], value=['not_expressed', 'expressed'])

# Merge predictions to feature df and save df
merged_df = feature_df_copy.merge(prediction_df, how='left', on='gene_id_sno')
merged_df.to_csv(snakemake.output.predicted_label_df, sep='\t', index=False)

# Save also scaled_feature_df and real label as dfs
scaled_feature_df = pd.DataFrame(scaled_feature_df, index=final_df.index, columns=final_df.columns)
scaled_feature_df = scaled_feature_df.reset_index()
scaled_feature_df.to_csv(snakemake.output.scaled_feature_df, sep='\t', index=False)
