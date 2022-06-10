#!/usr/bin/python3
import pandas as pd
import functions as ft
import collections as coll

df = pd.read_csv(snakemake.input.confusion_value_per_sno, sep='\t')
colors = snakemake.params.color_dict
output = snakemake.output.pie

col = list(df.filter(regex='^confusion_matrix_').columns)

# Generate a pie chart of the number of snoRNA per confusion_value (TP, TN, FP, FN)
counts = [len(df[df[col[0]] == 'TP']),
            len(df[df[col[0]] == 'TN']),
            len(df[df[col[0]] == 'FP']),
            len(df[df[col[0]] == 'FN'])]  # keep this order as in the config.json color_dict
ft.pie_simple(counts, colors, '', output)
