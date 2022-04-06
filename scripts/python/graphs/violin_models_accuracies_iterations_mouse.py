#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft

""" For each model (log_reg, svc and rf), represent a violin plot of the
    accuracies across the 10 iterations based on the predictions of the
    abundance status of mouse snoRNAs."""

accuracies_df_paths = snakemake.input.accuracies
color_dict = snakemake.params.color_dict

# Create dict containing the accuracy of all iterations per model 
acc = {}
for i, path in enumerate(accuracies_df_paths):
    model_name, iteration = path.split('/')[-1].split('.')[0].split('_test_accuracy_')
    df = pd.read_csv(path, sep='\t')
    accuracy = float(df.values[0][0])
    if model_name not in acc.keys():
        acc[model_name] = {iteration: accuracy}
    else:
        acc[model_name][iteration] = accuracy

# Create a df containing all the accuracies in the right format for the violin plot function
dfs = []
for mod_name, accuracy_dict in acc.items():
    vals = list(accuracy_dict.values())
    df = pd.DataFrame({mod_name: vals})
    df = df.rename(columns={mod_name: 'Accuracy'})
    df['Model'] = mod_name
    dfs.append(df)
    
final_df = pd.concat(dfs)

# Create the violin plot
ft.violin_wo_swarm(final_df, 'Model', 'Accuracy', None, 'Model', 'Accuracy of each\nmodel per iteration', '',
            color_dict, snakemake.output.violin)
