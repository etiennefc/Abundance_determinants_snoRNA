#!/usr/bin/python3
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)  # ignore all future warnings
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
import pickle

""" Train (fit) each model on the training set using the best
    hyperparameters found by hyperparameter_tuning_cv. Pickle these fitted
    models (into .sav files) so that they can be reused after without the
    need to retrain them all over again."""

# Get best hyperparameters per model
hyperparams_df = pd.read_csv(snakemake.input.best_hyperparameters, sep='\t')
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

# Generate the same CV, training and test sets (only the training set will be
# used in this script) that were generated in hyperparameter_tuning_cv
# (respectively 15%, 70% and 15% of all dataset examples)
df = pd.read_csv(snakemake.input.df, sep='\t', index_col='gene_id_sno')
X = df.drop('label', axis=1)
y = df['label']

# First the CV vs total_train split
X_total_train, X_cv, y_total_train, y_cv = train_test_split(X, y, test_size=0.15,
                                            random_state=42, stratify=y)

# Next the total_train is split into train and test sets (1017 and 180 correspond
# to the number of examples in train and test sets respectively to get an
# approximately 70 % and 15 % of all examples in these two datasets)
X_train, X_test, y_train, y_test = train_test_split(X_total_train, y_total_train,
                                    test_size=180, train_size=1017, random_state=42,
                                    stratify=y_total_train)


# Instantiate the model defined by the 'models' wildcard using the best hyperparameters
# specific to each model (log_reg, svc, rf, gbm, knn)
if snakemake.wildcards.models2 == "log_reg":
    model = LogisticRegression(C=hyperparams['C'], solver=hyperparams['solver'],
                                random_state=42, max_iter=500)
elif snakemake.wildcards.models2 == "svc":
    model = svm.SVC(C=hyperparams['C'], degree=hyperparams['degree'],
                    gamma=hyperparams['gamma'], kernel=hyperparams['kernel'],
                    random_state=42)
elif snakemake.wildcards.models2 == "rf":
    model = RandomForestClassifier(max_depth=hyperparams['max_depth'],
                min_samples_leaf=hyperparams['min_samples_leaf'],
                min_samples_split=hyperparams['min_samples_split'],
                n_estimators=hyperparams['n_estimators'], random_state=42)
elif snakemake.wildcards.models2 == "knn":
    model = KNeighborsClassifier(n_neighbors=hyperparams['n_neighbors'],
                weights=hyperparams['weights'],
                leaf_size=hyperparams['leaf_size'], p=hyperparams['p'])
else:
    model = GradientBoostingClassifier(loss=hyperparams['loss'],
                max_depth=hyperparams['max_depth'],
                min_samples_leaf=hyperparams['min_samples_leaf'],
                min_samples_split=hyperparams['min_samples_split'],
                n_estimators=hyperparams['n_estimators'], random_state=42)

# Train model and save training accuracy to df
model.fit(X_train, y_train)
print(snakemake.wildcards.models2)
print(model.score(X_train, y_train))

acc = {}
acc[snakemake.wildcards.models2+'_training_accuracy'] = model.score(X_train, y_train)
acc_df = pd.DataFrame(acc, index=[0])
acc_df.to_csv(snakemake.output.training_accuracy, sep='\t', index=False)

# Pickle the model as a .sav file ('wb' for write in binary)
pickle.dump(model, open(snakemake.output.pickled_trained_model, 'wb'))
