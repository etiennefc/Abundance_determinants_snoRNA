#!/usr/bin/python3
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.utils import shuffle

""" Fill NaN in feature df in numerical columns with -5 instead. -5 was chosen
    arbitrarily so that this negative value should not interfere with all other
    positive values. Then, AFTER splitting into train, CV and test sets, apply
    mean normalization (feature scaling) to numerical features columns."""
all_iterations = ["manual_first", "manual_second", "manual_third", "manual_fourth",
                "manual_fifth", "manual_sixth", "manual_seventh", "manual_eighth",
                "manual_ninth", "manual_tenth"]
manual_iteration = snakemake.wildcards.manual_iteration
idx = all_iterations.index(manual_iteration)
random_state = snakemake.params.random_state
df_0 = pd.read_csv(snakemake.input.feature_df, sep='\t')

# Fill NaN with -5
df_0 = df_0.fillna(-5)


# First, shuffle the df
df = df_0.copy()
df = shuffle(df, random_state=random_state)

X = df[['gene_id_sno', 'combined_box_hamming']]
y = df['label']

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
        y_test.to_csv(snakemake.output.y_test, index=False, sep='\t')
        y_train.to_csv(snakemake.output.y_train, index=False, sep='\t')
        y_cv.to_csv(snakemake.output.y_cv, index=False, sep='\t')

        # Scale feature values using mean normalization for numerical value columns
        # with high standard deviation
        dfs = [X_cv, X_train, X_test]
        output = [snakemake.output.cv, snakemake.output.train, snakemake.output.test]
        for j, df in enumerate(dfs):
            df_num = df.select_dtypes(include=['int64', 'float64'])
            num_cols = list(df_num.columns)
            for i, col in enumerate(num_cols):
                mean = df[col].mean()
                std = df[col].std()
                if std != 0:
                    df[col+'_norm'] = (df[col] - mean) / std
                else:  # to deal with column that has all the same value, thus a std=0
                    df[col+'_norm'] = df[col]  # we don't scale, but these values will either be all 0 or all 1
            df = df.drop(num_cols, axis=1)
            df.to_csv(output[j], index=False, sep='\t')
    i+=1
