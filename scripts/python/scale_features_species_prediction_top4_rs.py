#!/usr/bin/python3
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler

""" Fill NaN in feature df in numerical columns with -5 instead. -5 was chosen
    arbitrarily so that this negative value should not interfere with all other
    positive values. Then, AFTER splitting into train and CV sets, apply
    mean normalization (feature scaling) to numerical features columns."""

random_state = int(snakemake.wildcards.rs)

df_0 = pd.read_csv(snakemake.input.feature_df, sep='\t')

# Fill NaN with -5
df_0 = df_0.fillna(-5)


# First, shuffle the df
df = df_0.copy()
df = shuffle(df, random_state=random_state)

# Keep only top 4 features
X = df.drop('label', axis=1)
X = X[['gene_id_sno', 'combined_box_hamming', 'sno_mfe', 'terminal_stem_mfe', 'host_expressed']]
y = df['label']


# Split the X dataset into cv and train test (respectively 10 and 90 % of the total df)
X_train, X_cv, y_train, y_cv = train_test_split(X, y,
                            test_size=0.1, train_size=0.9, random_state=random_state,
                            stratify=y)
y_cv.to_csv(snakemake.output.y_cv, index=False, sep='\t')
y_train.to_csv(snakemake.output.y_train, index=False, sep='\t')

# Scale feature values using mean normalization for numerical value columns
dfs = [X_cv.set_index('gene_id_sno'), X_train.set_index('gene_id_sno')]
output = [snakemake.output.cv, snakemake.output.train]
for j, df in enumerate(dfs):
    scaler = StandardScaler().fit(df)
    scaled_df_array = scaler.transform(df)
    scaled_df = pd.DataFrame(scaled_df_array, index=df.index, columns=df.columns)
    scaled_df = scaled_df.reset_index()
    scaled_df.to_csv(output[j], index=False, sep='\t')
