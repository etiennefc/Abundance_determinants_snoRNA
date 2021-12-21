#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import functions as ft
import statistics as sts

""" Generate an upset plot per confusion value (true positive, false positive,
    false negative or true negative) to see the intersection in the snoRNAs
    (mis)classified by all three models not overfitting (log_reg, svc and rf).
    The upset plot shows the average intersection (+/-stdev) across 10 iterations
    for all models."""
df_output_path = snakemake.params.df_output_path

# Sort alphabetically all confusion matrix df (one per iteration) per model
log_reg_conf = snakemake.input.log_reg
log_reg_conf = sorted(log_reg_conf)
svc_conf = snakemake.input.svc
svc_conf = sorted(svc_conf)
rf_conf = snakemake.input.rf
rf_conf = sorted(rf_conf)

# Load confusion matrix of all models per iteration
log_reg_dfs, svc_dfs, rf_dfs = [], [], []
for i, log_reg_df in enumerate(log_reg_conf):
    log_reg = pd.read_csv(log_reg_conf[i], sep='\t')
    log_reg = log_reg[['gene_id_sno', 'confusion_matrix_val_log_reg']]
    log_reg_dfs.append(log_reg)
    svc = pd.read_csv(svc_conf[i], sep='\t')
    svc = svc[['gene_id_sno', 'confusion_matrix_val_svc']]
    svc_dfs.append(svc)
    rf = pd.read_csv(rf_conf[i], sep='\t')
    rf = rf[['gene_id_sno', 'confusion_matrix_val_rf']]
    rf_dfs.append(rf)

# Get the avg (and stdev) total number of confusion_value (ex: TN) across iterations per model
# This will be used for the horizontal bar chart within the upset plot
confusion_val_nb_per_iteration_log_reg = []
confusion_val_nb_per_iteration_svc = []
confusion_val_nb_per_iteration_rf = []
for i, log_reg_df in enumerate(log_reg_dfs):
    nb_log_reg = len(log_reg_df[log_reg_df['confusion_matrix_val_log_reg'] == snakemake.wildcards.confusion_value])
    confusion_val_nb_per_iteration_log_reg.append(nb_log_reg)
    nb_svc = len(svc_dfs[i][svc_dfs[i]['confusion_matrix_val_svc'] == snakemake.wildcards.confusion_value])
    confusion_val_nb_per_iteration_svc.append(nb_svc)
    nb_rf = len(rf_dfs[i][rf_dfs[i]['confusion_matrix_val_rf'] == snakemake.wildcards.confusion_value])
    confusion_val_nb_per_iteration_rf.append(nb_rf)

# Order from rf to log_reg (bottom to top of horizontal bar chart)
avg_nb, stdev_nb = [], []
for model in [confusion_val_nb_per_iteration_rf, confusion_val_nb_per_iteration_svc, confusion_val_nb_per_iteration_log_reg]:
    avg_, stdev_ = sts.mean(model), sts.stdev(model)
    avg_nb.append(avg_)
    stdev_nb.append(stdev_)


# This will be used for the scatter and vertical bar plot in the upset plot
# Merge all dfs within each iteration to create one df
# and get a dict per iteration to show the number of elements per upset plot
# category (e.g. TP_TP_FN: 32, (TP for log_reg, then TP for svc and finally FN for rf))
upset_dicts = []
temp_dfs = []
for i, log_reg_df in enumerate(log_reg_dfs):
    df = log_reg_dfs[i].merge(svc_dfs[i], how='left', on='gene_id_sno')
    df = df.merge(rf_dfs[i], how='left', on='gene_id_sno')
    df.to_csv(df_output_path[i], sep='\t', index=False)


    # Select rows containing at least one confusion_value wildcard (ex: TP)
    val_df = df[(df.iloc[:, 1:4] == snakemake.wildcards.confusion_value).any(axis=1)]
    unique_category = val_df.drop('gene_id_sno', axis=1)
    unique_category = unique_category.drop_duplicates(['confusion_matrix_val_log_reg',
                            'confusion_matrix_val_svc',
                            'confusion_matrix_val_rf'])
    temp_dfs.append(unique_category)  # get only unique combination of confusion value per iteration (merged_df of that iteration)

    # Get the occurences of all possible categories containing the confusion_value wildcard (ex: TP_TP_FN, TP_TP_TP, etc.)
    groups = val_df.groupby(['confusion_matrix_val_log_reg',
                            'confusion_matrix_val_svc',
                            'confusion_matrix_val_rf'])
    d = {}
    for i, group in enumerate(groups):
        upset_category_name = "_".join(group[0])
        number_per_category = len(group[1])
        d[upset_category_name] = number_per_category
    upset_dicts.append(d)

# Get the union of all upset possible categories (ex: TP_TP_FN, TP_TP_TP, etc.) across iterations
concat_temp_df = pd.concat(temp_dfs)
concat_temp_df['upset_category_name'] = concat_temp_df['confusion_matrix_val_log_reg'] + '_' + concat_temp_df['confusion_matrix_val_svc'] + '_' + concat_temp_df['confusion_matrix_val_rf']
union = list(pd.unique(concat_temp_df['upset_category_name']))

# Add 0 as value to missing keys in each dictionary in the upset_dicts list (so that all iterations have the same number of categories)
for cat in union:
    for d in upset_dicts:
        d.setdefault(cat, 0)

# Get the average and stdev for all upset categories across iterations
upset_df = pd.DataFrame(upset_dicts)
average = {}
stdev = {}
for col in upset_df.columns:
    avg, std = upset_df[col].mean(), upset_df[col].std()
    average[col], stdev[col] = avg, std

# Sort by descending order of value
sorted_average = sorted(average.items(), key=lambda x: x[1], reverse=True)
sorted_stdev = sorted(stdev.items(), key=lambda x: x[1], reverse=True)


# Convert confusion value names to model names (ex: if TP for TP_FN_TP, then it
# will output the first and third model name, but not the second) for scatter plot in the upset
def convert_vals(val_list, confusion_value):
    """ ex. of val_list --> ['TP_TP_TP', 'TP_TN_TN', etc.]
        ex. of confusion_value --> 'TP'"""
    name, values = [], []
    for val in val_list:
        log_reg, svc, rf = val.split('_')
        for i, confusion_val in enumerate([rf, svc, log_reg]):  # the order is from y=0 to y=2 on the scatter
            if confusion_val == confusion_value:
                name.append(val)
                values.append(i)  # dot at y=i on the scatter
    return name, values

names, values = convert_vals([key[0] for key in sorted_average], snakemake.wildcards.confusion_value)


# Get minimal and maximal values (0, 1 or 2) of each upset category and their position on the x axis
ymins, ymaxs, vlines_pos = [], [], []
for name in list(set(names)):
    index_in_names = [i for i in range(len(names)) if names[i] == name]
    vals = [values[j] for j in index_in_names]

    # get position of each vertical line on the x-axis (according to sorted_average, i.e values for the bar chart sorted in descending order)
    pos = [sorted_average.index(key) for key in sorted_average if key[0] == name] # list of only one element

    # get minimal and maximal values for each vertical line on the upset plot
    mini, maxi = min(vals), max(vals)
    ymins.append(mini)
    ymaxs.append(maxi)
    vlines_pos.append(pos[0])

# Create homemade upset plot (vertical bar chart over dot plot plus a horizontal bar chart)
ft.upset_avg_3_cat(sorted_average, sorted_stdev, avg_nb, stdev_nb, names, values,
                    vlines_pos, ymins, ymaxs, 'Average intersection size\nacross iterations',
                    'Average number\nper model\nacross iterations',
                    ["RandomForest", "SupportVector", "LogisticRegression"],
                    snakemake.wildcards.confusion_value, snakemake.output.upset)
