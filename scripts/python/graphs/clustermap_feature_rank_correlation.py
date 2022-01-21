#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft
import numpy as np
import seaborn as sns

rank_df = pd.read_csv(snakemake.input.rank_features_df, sep='\t')
rank_df[['feature2', 'norm']] = rank_df['feature'].str.split('_norm', expand=True)
rank_df = rank_df.drop(columns=['norm', 'feature'])

feature_distribution = {}
for i, group in enumerate(rank_df.groupby('feature2')['feature_rank']):
    feature_name = group[0]
    range_ = group[1].max() - group[1].min()
    median_ = group[1].median()
    feature_distribution[feature_name] = [median_, range_]

# Order features by increasing median value of feature_ranks and by range as second sort if two features have the same median
# i.e. the same order as in the violin plot of feature rank
feature_distribution_df = pd.DataFrame.from_dict(feature_distribution, columns = ['median', 'range'], orient='index')
ordered_features = feature_distribution_df.sort_values(by=['median', 'range'], ascending=[True, True]).index.to_list()
print(ordered_features)


models = ['log_reg', 'svc', 'rf']
outputs = snakemake.output.heatmaps
for mod in models:
    output = [path for path in outputs if mod in path][0]
    iterations_df = rank_df[rank_df['model'].str.startswith(mod)]
    print(list(iterations_df['feature2']))
    pivot = iterations_df.pivot(index='model', columns='feature2', values='feature_rank')
    pivot = pivot[ordered_features]
    print(pivot)
    correlation_df = pivot.corr(method='spearman')
    print(correlation_df)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    sns.clustermap(correlation_df, cmap='viridis', cbar_kws={'label': "Feature rank correlation\n(Spearman's ρ)"},
                row_cluster=False)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.xlabel(xlabel="Features")
    plt.ylabel(ylabel="Features")
    plt.savefig(output, dpi=600, bbox_inches='tight')
