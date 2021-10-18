#!/usr/bin/python3
import pandas as pd
from pandas.core.common import flatten
import functions as ft
import collections as coll

df = pd.read_csv(snakemake.input.df, sep='\t')
colors = snakemake.params.rank_colors
output = snakemake.output.donut

def count_in_list(input_list, ref_list):
    """ Count the number of occurences of each element in ref_list within the
        input_list. Return the associated percentages in a list."""
    l = []
    for element in ref_list:
        occurences = input_list.count(element)
        l.append(occurences)

    total = sum(l)
    percent_list = []
    for count in l:
        if total != 0:
            percent = count/total * 100
            percent = round(percent, 2)
            percent_list.append(percent)
        else:
            percent = 0
            percent_list.append(percent)

    return percent_list


# Create a dict where each key is a feature and each value is a list of all ranks
# (corresponding to that feature) within the top 5 predictive features across models that use that feature
rank_dict = df.groupby('feature')['feature_rank'].apply(list).to_dict()
print(rank_dict)

## Create a dict where the keys are the models intersections
## (e.g., gbm, log_reg_svc_gbm_knn, svc_knn, etc.) and the values are the features shared by these models
# First aggregate model names together (ex: if a feature was shared by svc and knn, then the model name becomes svc_knn)
groups = df.groupby('feature').agg({'model':list})
groups['model'] = groups['model'].apply(lambda x: "_".join(map(str, x)))
groups.reset_index(inplace=True)
# Combine the common features to each models intersection into a list and create the dict
groups = groups.groupby('model')['feature'].apply(lambda x: ",".join(map(str, x))).reset_index()
groups['feature'] = groups['feature'].str.split(',')
feature_model_intersect_dict = dict(zip(groups.model, groups.feature))

print(feature_model_intersect_dict)

# Generate list of list of rank percentage (1 list per intersection category (8 categories in total))
model_intersect_names = []
features_name = []
all_rank_percent = []
for model_intersect, features in feature_model_intersect_dict.items():
    if len(features) == 1:  # if only one feature is shared in that model intersection
        ranks = rank_dict[features[0]]
        rank_percent = count_in_list(ranks, [1,2,3,4,5])  # count how many rank 1,2, ..., 5 are associated to that feature across model(s)
        all_rank_percent.append(rank_percent)
    elif len(features) > 1:  # if multiple features are common in that model intersection
        temp = []
        for feature in features:
            ranks = rank_dict[feature]
            temp.append(ranks)
        rank_percent = count_in_list(list(flatten(temp)), [1,2,3,4,5])
        all_rank_percent.append(rank_percent)
    model_intersect_names.append(model_intersect)
    features_name.append(features)

print(model_intersect_names)
print(features_name)
print(all_rank_percent)

# Generate a donut chart per model intersection according to their rank (1st, 2nd, ..., 5th most predictive feature) across models
ft.pie_multiple(all_rank_percent, colors.keys(), colors.values(), model_intersect_names,
    'Ranking of top 5 most predictive features across model intersections', 'Rank across top 5 features and models', output)
