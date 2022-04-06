#!/usr/bin/python3
import pandas as pd
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.utils import shuffle

""" Scale features based on the training sets used to train the models and
    predict the abundance status of SNORA81s."""
all_iterations = ["manual_first", "manual_second", "manual_third", "manual_fourth",
                "manual_fifth", "manual_sixth", "manual_seventh", "manual_eighth",
                "manual_ninth", "manual_tenth"]
manual_iteration = snakemake.wildcards.manual_iteration
idx = all_iterations.index(manual_iteration)
random_state = snakemake.params.random_state
human_feature_df_0 = pd.read_csv(snakemake.input.human_snoRNA_feature_df, sep='\t')
feature_df = pd.read_csv(snakemake.input.feature_df, sep='\t', index_col='gene_id_sno')

# Unpickle and thus instantiate the trained model defined by the 'models' wildcard
model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))

# Fill NaN with -5
human_feature_df_0 = human_feature_df_0.fillna(-5)
human_feature_df = human_feature_df_0.copy()
human_feature_df = shuffle(human_feature_df, random_state=random_state)

# Keep only top 3 features
X = human_feature_df.drop('label', axis=1)
X = X[['gene_id_sno', 'sno_mfe', 'terminal_stem_mfe', 'combined_box_hamming', 'host_expressed']]  # same order as in feature_df that we want to predict
y = human_feature_df['label']

# Configure the cross-validation strategy (StratifiedKFold where k=10)
# This serves only to split in ten equal folds, no cross-validation is done to this point
kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=random_state)

# Get the corresponding test set (10% of total df) for that manual iteration
i = 0
for total_train_index, test_index in kf.split(X, y):
    if i == idx:
        X_test = X.loc[test_index]
        y_test = y.loc[test_index]

        # Split the total_train into cv and train test (respectively 10 and 80 % of the total df (i.e. 11.1% and 88.9% of 1386 snoRNAs))
        X_total_train = X.loc[total_train_index]
        y_total_train = y.loc[total_train_index]
        X_train, X_cv, y_train, y_cv = train_test_split(X_total_train, y_total_train,
                                    test_size=0.111, train_size=0.889, random_state=random_state,
                                    stratify=y_total_train)
        X_train = X_train.drop('gene_id_sno', axis=1)
        # Scale by substracting mean and dividing by stdev
        scaler = StandardScaler().fit(X_train)
    i+=1

# Scale new feature_df (SNORA81s) by the same scaling factor used to scale the training set
scaled_feature_df = scaler.transform(feature_df)

# Predict label (expressed (1) or not_expressed (0)) on test data and compare to y_test
y_pred = model.predict(scaled_feature_df)
prediction_df = pd.DataFrame(y_pred, index=feature_df.index, columns=['label'])
prediction_df = prediction_df.reset_index()
prediction_df.to_csv(snakemake.output.predicted_label_df, sep='\t', index=False)
