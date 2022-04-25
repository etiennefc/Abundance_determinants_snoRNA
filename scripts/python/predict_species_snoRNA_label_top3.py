#!/usr/bin/python3
import pandas as pd
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.utils import shuffle
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder

""" Scale features based on the training sets used to train the models and
    predict the abundance status of mouse snoRNAs."""

random_state = snakemake.params.random_state
human_feature_df_0 = pd.read_csv(snakemake.input.human_snoRNA_feature_df, sep='\t')
feature_df = pd.read_csv(snakemake.input.feature_df, sep='\t')  # mouse snoRNAs


## Mouse dataset proocessing
# Label-encode the label column manually
feature_df.loc[feature_df['abundance_cutoff'] == 'expressed', 'label'] = 1
feature_df.loc[feature_df['abundance_cutoff'] == 'not_expressed', 'label'] = 0
feature_df = feature_df.drop(columns=['gene_name', 'abundance_cutoff'])

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
y_mouse = final_df['label']
final_df = final_df[['sno_mfe', 'terminal_stem_mfe', 'combined_box_hamming']]


## Human dataset processing
# Fill NaN with -5
human_feature_df_0 = human_feature_df_0.fillna(-5)
human_feature_df = human_feature_df_0.copy()
human_feature_df = shuffle(human_feature_df, random_state=random_state)

# Keep only top 3 features
X = human_feature_df.drop('label', axis=1)
X = X[['gene_id_sno', 'sno_mfe', 'terminal_stem_mfe', 'combined_box_hamming']]  # same order as in final_df that we want to predict
y = human_feature_df['label']


# Split the X dataset into cv and train test (respectively 10 and 90 % of the total df)
X_train, X_cv, y_train, y_cv = train_test_split(X, y,
                            test_size=0.1, train_size=0.9, random_state=random_state,
                            stratify=y)

# Create the scaler that will be used to mean normalize mouse snoRNA features
df = X_train.set_index('gene_id_sno')
scaler = StandardScaler().fit(df)



# Scale final_df (all mouse snoRNAs with top3 features) by the same scaling factor used to scale the training set in human snoRNAs
scaled_feature_df = scaler.transform(final_df)

# Unpickle and thus instantiate the trained model defined by the 'models' wildcard
model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))

# Predict label (expressed (1) or not_expressed (0)) on mouse snoRNAs
y_pred = model.predict(scaled_feature_df)
prediction_df = pd.DataFrame(y_pred, index=final_df.index, columns=['predicted_label'])
prediction_df = prediction_df.reset_index()
prediction_df.to_csv(snakemake.output.predicted_label_df, sep='\t', index=False)

# Save also scaled_feature_df and real label as dfs
scaled_feature_df = pd.DataFrame(scaled_feature_df, index=final_df.index, columns=final_df.columns)
scaled_feature_df = scaled_feature_df.reset_index()
scaled_feature_df.to_csv(snakemake.output.scaled_feature_df, sep='\t', index=False)
y_mouse.to_csv(snakemake.output.label_df, sep='\t', index=False)
