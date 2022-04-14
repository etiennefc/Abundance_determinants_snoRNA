#!/usr/bin/python3
import pandas as pd
import collections as coll
from statistics import mode

""" Define the consensus confusion value across the 10 iterations for each model.
    Choose the confusion value based on the highest number of time predicted
    as such across the 10 models. If 2 confusion values have an equal number of
    votes (5 vs 5), remove randomly 1 iteration and the equality will then be
    broken and a predominant confusion value will be chosen."""

df_paths = snakemake.input.confusion_val_df

all_mouse_sno_ids = pd.read_csv(df_paths[0], sep='\t')  # we take the first df here, but they all contain all the snoRNAs
all_mouse_sno_ids = list(all_mouse_sno_ids['gene_id_sno'])

conf_val = {}
for sno_id in all_mouse_sno_ids:
    temp_val = []
    for path in df_paths:
        df = pd.read_csv(path, sep='\t')
        df = df.filter(regex='gene_id_sno|^confusion_matrix_val')
        val = df[df['gene_id_sno'] == sno_id].values[0][1]
        temp_val.append(val)
    if 5 in coll.Counter(temp_val).values():  # if there is an equality in confusion value votes (ex: 5 TN vs 5 FP)
        del temp_val[0]  # remove first iteration 
        consensus = mode(temp_val)
        conf_val[sno_id] = consensus
    else:
        consensus = mode(temp_val)
        conf_val[sno_id] = consensus

# Create df from dict
final_df = pd.DataFrame(conf_val.items(), columns=['gene_id_sno', 'consensus_confusion_value'])
final_df.to_csv(snakemake.output.consensus_conf_val_df, sep='\t', index=False)
