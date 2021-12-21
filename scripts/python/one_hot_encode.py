#!/usr/bin/python3
import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder

""" One-hot encode categorical features and label-encode the labels."""

df = pd.read_csv(snakemake.input.feature_df, sep='\t')

# Label-encode the label column manually
df.loc[df['abundance_cutoff_2'] == 'expressed', 'label'] = 1
df.loc[df['abundance_cutoff_2'] == 'not_expressed', 'label'] = 0
df = df.drop(columns=['gene_name', 'abundance_cutoff', 'abundance_cutoff_2'])

numerical_features_label = df.copy()
numerical_features_label = numerical_features_label.select_dtypes(include=['int64', 'float64'])

# Convert categorical features (astype 'object') into numpy array, convert the
# strings in that array into numbers with LabelEncoder and then one-hot encode this numerical array
dfs = [df[['gene_id_sno']]]
df = df.drop(columns=['gene_id_sno'])
df_cat = df.select_dtypes(include=['object'])
cols = df_cat.columns
for i, col in enumerate(cols):
    df_cat = df[[col]]

    # Convert column into numpy array
    array_cat = df_cat.values.reshape(-1, 1)  # -1 infers the length of df_cat (i.e. 1541)

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

# Remove duplicated columns to keep only one (e.g. 'intergenic', which is created 5 times when one-hot encoding host-related columns)
final_df = final_df.loc[:,~final_df.columns.duplicated()]

final_df.to_csv(snakemake.output.one_hot_encoded_df, sep='\t')
