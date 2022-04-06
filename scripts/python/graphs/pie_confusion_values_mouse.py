#!/usr/bin/python3
import pandas as pd
import functions as ft

df = pd.read_csv(snakemake.input.confusion_value_per_sno, sep='\t')
colors = snakemake.params.color_dict
output = snakemake.output.pie

# Generate a pie chart of the number of snoRNA per confusion_value (TP, TN, FP, FN)
counts = [len(df[df['consensus_confusion_value'] == 'TP']),
            len(df[df['consensus_confusion_value'] == 'TN']),
            len(df[df['consensus_confusion_value'] == 'FP']),
            len(df[df['consensus_confusion_value'] == 'FN'])]  # keep this order as in the config.json color_dict
ft.pie_simple(counts, colors, '', output)
