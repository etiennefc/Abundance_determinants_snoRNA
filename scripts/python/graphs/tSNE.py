#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import numpy as np

"""
Creates a tSNE plot (t-distributed Stochastic Neighbor Embedding).
"""

df_initial = pd.read_csv(snakemake.input.df, sep='\t')
y = df_initial['label']

# First the CV vs total_train split
X_total_train, X_cv, y_total_train, y_cv = train_test_split(df_initial, y, test_size=0.15,
                                            random_state=42, stratify=y)

# Next the total_train is split into train and test sets (1077 and 232 correspond
# to the number of examples in train and test sets respectively to get an
# approximately 70 % and 15 % of all examples in these two datasets)
X_train, X_test, y_train, y_test = train_test_split(X_total_train, y_total_train,
                                    test_size=232, train_size=1077, random_state=42,
                                    stratify=y_total_train)


X = df_initial.drop(['label', 'gene_id_sno'], axis=1)
X_test_copy = X_test.copy()
X_test_copy = X_test_copy.drop(['label', 'gene_id_sno'], axis=1)


# Normalize data for tSNE
val = X.values
norm_val = StandardScaler().fit_transform(val)

val_test = X_test_copy.values
norm_val_test = StandardScaler().fit_transform(val_test)


# Create tSNE with 2 components for all snoRNAs or only those in the test set
all_sno_t_sne = TSNE(n_components=2, random_state=42).fit_transform(norm_val)
df_initial['Component_1'] = all_sno_t_sne[:, 0]
df_initial['Component_2'] = all_sno_t_sne[:, 1]
df_initial['label'] = df_initial['label'].replace([0, 1], ['not_expressed', 'expressed'])
df_initial['intergenic'] = df_initial['intergenic'].replace([0, 1], ['intronic', 'intergenic'])

test_sno_t_sne = TSNE(n_components=2, random_state=42).fit_transform(norm_val_test)
X_test['Component_1'] = test_sno_t_sne[:, 0]
X_test['Component_2'] = test_sno_t_sne[:, 1]
X_test['label'] = X_test['label'].replace([0, 1], ['not_expressed', 'expressed'])
X_test['intergenic'] = X_test['intergenic'].replace([0, 1], ['intronic', 'intergenic'])

# Create the (tSNE) scatter plot for all snoRNAs
plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(15, 15))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax.set_xlabel('Dimension 1', fontsize=20)
ax.set_ylabel('Dimension 2', fontsize=20)

colors = list(snakemake.params.colors_dict.values())
sns.scatterplot(x='Component_1', y='Component_2', data=df_initial,
                hue=snakemake.wildcards.pca_hue, palette=colors, ax=ax)
plt.savefig(snakemake.output.t_sne_all, dpi=600)


# Create the (tSNE) scatter plot for snoRNAs in test set only
plt.rcParams['svg.fonttype'] = 'none'
fig2, ax2 = plt.subplots(1, 1, figsize=(15, 15))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax2.set_xlabel('Dimension 1', fontsize=20)
ax2.set_ylabel('Dimension 2', fontsize=20)

colors = list(snakemake.params.colors_dict.values())
sns.scatterplot(x='Component_1', y='Component_2', data=X_test,
                hue=snakemake.wildcards.pca_hue, palette=colors, ax=ax2)
plt.savefig(snakemake.output.t_sne_test, dpi=600)
