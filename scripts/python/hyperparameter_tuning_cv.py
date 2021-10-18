#!/usr/bin/python3
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV, RandomizedSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

""" Tune the hyperparameters of each models (Logistic regression, Support
    vector classifier, Random Forest and Gradient boosting classifier)
    before even training them, using GridSearchCV with stratified k-fold.
    First, shuffle the dataset and split it into Cross-validation (CV) set
    and training set (respectively 15 % and 85 % of all examples). The CV
    will take place as a stratifed k-fold (3 fold) using GridSearchCV and
    will return the best hyperparameters for each tested model. Of note, the
    training set will be also split later on into training and test set
    (respectively 70 % and 15 % of all examples)."""

df = pd.read_csv(snakemake.input.df, sep='\t', index_col='gene_id_sno')
X = df.drop('label', axis=1)
y = df['label']

# Split dataset into CV set (15 % of all examples) and training set (85 %) and use
# the 'stratify' param to keep the same proportion of expressed vs not_expressed in
# training and CV sets. The datasets are shuffled by default
X_train, X_cv, y_train, y_cv = train_test_split(X, y, test_size=0.15, random_state=42, stratify=y)

# Instantiate the model defined by the 'models' wildcard
if snakemake.wildcards.models2 == "log_reg":
    model = LogisticRegression(random_state=42, max_iter=500)
elif snakemake.wildcards.models2 == "svc":
    model = svm.SVC(random_state=42)
elif snakemake.wildcards.models2 == "rf":
    model = RandomForestClassifier(random_state=42)
elif snakemake.wildcards.models2 == "knn":
    model = KNeighborsClassifier()
else:
    model = GradientBoostingClassifier(random_state=42)

# Configure the cross-validation strategy (StratifiedKFold where k=3)
cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

# Define search space (i.e. in which range of params the GridSearch happens) per model
space = snakemake.params.hyperparameter_space

# Execute the gridsearch per model on the CV set
search = GridSearchCV(estimator=model, param_grid=space,
                        cv=cv, scoring="accuracy")
search.fit(X_cv, y_cv)
print(snakemake.wildcards.models2)
print(search.best_score_)
print(search.best_params_)

# Return the best hyperparameters found by GridSearchCV, and the accuracy of each model
# fitted on the CV set with these hyperparameters into a dataframe
params_dict = search.best_params_
params_dict['accuracy_cv'] = search.best_score_
params_df = pd.DataFrame(params_dict, index=[0])
params_df.to_csv(snakemake.output.best_hyperparameters, sep='\t', index=False)
