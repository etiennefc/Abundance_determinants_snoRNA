#!/usr/bin/python3
import pandas as pd
from sklearn.model_selection import train_test_split

""" Fill NaN in feature df in numerical columns with -5 instead. -5 was chosen
    arbitrarily so that this negative value should not interfere with all other
    positive values. Then, AFTER splitting into train, CV and test sets, apply
    mean normalization (feature scaling) to numerical features columns."""
iteration = snakemake.wildcards.iteration_20
random_state_dict = snakemake.params.random_state
random_state = random_state_dict[iteration]
df = pd.read_csv(snakemake.input.feature_df, sep='\t')

# Fill NaN with -5
df = df.fillna(-5)

X = df.drop('label', axis=1)
y = df['label']

# First the CV vs total_train split
X_total_train, X_cv, y_total_train, y_cv = train_test_split(X, y, test_size=0.15,
                                            random_state=random_state, stratify=y)
y_cv.to_csv(snakemake.output.y_cv, index=False, sep='\t')

# Next the total_train is split into train and test sets (1077 and 232 correspond
# to the number of examples in train and test sets respectively to get an
# approximately 70 % and 15 % of all examples in these two datasets)
X_train, X_test, y_train, y_test = train_test_split(X_total_train, y_total_train,
                                    test_size=232, train_size=1077, random_state=random_state,
                                    stratify=y_total_train)
y_train.to_csv(snakemake.output.y_train, index=False, sep='\t')
y_test.to_csv(snakemake.output.y_test, index=False, sep='\t')

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
        df[col+'_norm'] = (df[col] - mean) / std
    df = df.drop(num_cols, axis=1)
    df.to_csv(output[j], index=False, sep='\t')
