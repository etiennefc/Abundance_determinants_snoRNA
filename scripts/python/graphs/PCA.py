#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import decomposition as dec
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import numpy as np

"""
Creates a PCA plot (principal component analysis).
"""

df_initial = pd.read_csv(snakemake.input.df, sep='\t')
df_copy = df_initial.copy()
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


# Normalize data for PCA
val = X.values
norm_val = StandardScaler().fit_transform(val)

val_test = X_test_copy.values
norm_val_test = StandardScaler().fit_transform(val_test)


# Create PCA analysis with 2 principal components for all snoRNAs
pca_all = dec.PCA(n_components=2, random_state=42)
principal_components_all = pca_all.fit_transform(norm_val)
principal_df_all = pd.DataFrame(data=principal_components_all, columns=['Principal_component_1', 'Principal_component_2'])
print('Explained variation per principal component: {}'.format(pca_all.explained_variance_ratio_))  # returns the proportion of variance explained by each component
print('For each component, the proportion of each columns composing the component is: {}'.format(pca_all.components_))  # returns an array of the proportion that each column contributes per component;

# Create PCA analysis with 2 principal components for snoRNAs in test set only
pca_test = dec.PCA(n_components=2, random_state=42)
principal_components_test = pca_test.fit_transform(norm_val_test)
principal_df_test = pd.DataFrame(data=principal_components_test, columns=['Principal_component_1', 'Principal_component_2'])
print('Explained variation per principal component: {}'.format(pca_test.explained_variance_ratio_))
print('For each component, the proportion of each columns composing the component is: {}'.format(pca_test.components_))


# Create the (pca) scatter plot for all snoRNAs
pc1_all = round(pca_all.explained_variance_ratio_[0] * 100, 2)
pc2_all = round(pca_all.explained_variance_ratio_[1] * 100, 2)

plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(15, 15))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax.set_xlabel(f'Principal Component 1 ({pc1_all} %)', fontsize=15)
ax.set_ylabel(f'Principal Component 2 ({pc2_all} %)', fontsize=15)
principal_df2_all = pd.concat([principal_df_all, df_copy[['label', 'intergenic']]], axis=1)
principal_df2_all['label'] = principal_df2_all['label'].replace([0, 1], ['not_expressed', 'expressed'])
principal_df2_all['intergenic'] = principal_df2_all['intergenic'].replace([0, 1], ['intronic', 'intergenic'])

crits = list(snakemake.params.colors_dict.keys())
colors = list(snakemake.params.colors_dict.values())

# Plot each hue separately on the same ax
for crit, color in zip(crits, colors):
    indicesToKeep = principal_df2_all[snakemake.wildcards.pca_hue] == crit
    ax.scatter(principal_df2_all.loc[indicesToKeep, 'Principal_component_1'],
                principal_df2_all.loc[indicesToKeep, 'Principal_component_2'],
                c=color, s=50)

plt.legend(crits, prop={'size': 15})
plt.savefig(snakemake.output.pca_all, dpi=600)


# Create the (pca) scatter plot for snoRNAs in test set only
pc1_test = round(pca_test.explained_variance_ratio_[0] * 100, 2)
pc2_test = round(pca_test.explained_variance_ratio_[1] * 100, 2)

plt.rcParams['svg.fonttype'] = 'none'
fig2, ax2 = plt.subplots(1, 1, figsize=(15, 15))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax2.set_xlabel(f'Principal Component 1 ({pc1_test} %)', fontsize=15)
ax2.set_ylabel(f'Principal Component 2 ({pc2_test} %)', fontsize=15)
principal_df2_test = pd.concat([principal_df_test, X_test[['label', 'intergenic']].reset_index()], axis=1)
principal_df2_test['label'] = principal_df2_test['label'].replace([0, 1], ['not_expressed', 'expressed'])
principal_df2_test['intergenic'] = principal_df2_test['intergenic'].replace([0, 1], ['intronic', 'intergenic'])

crits = list(snakemake.params.colors_dict.keys())
colors = list(snakemake.params.colors_dict.values())


# Plot each hue separately on the same ax
for crit, color in zip(crits, colors):
    indicesToKeep = principal_df2_test[snakemake.wildcards.pca_hue] == crit
    ax2.scatter(principal_df2_test.loc[indicesToKeep, 'Principal_component_1'],
                principal_df2_test.loc[indicesToKeep, 'Principal_component_2'],
                c=color, s=50)

plt.legend(crits, prop={'size': 15})
plt.savefig(snakemake.output.pca_test, dpi=600)
