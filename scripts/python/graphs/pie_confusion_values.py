#!/usr/bin/python3
import pandas as pd
import functions as ft
import numpy as np

df = pd.read_csv(snakemake.input.sno_per_confusion_value, sep='\t')
colors = snakemake.params.color_dict
output = snakemake.output.pie
# Generate a pie chart of the number of snoRNA per confusion_value (TP, TN, FP, FN)
counts = [len(df[df['confusion_matrix'] == 'TP']),
            len(df[df['confusion_matrix'] == 'TN']),
            len(df[df['confusion_matrix'] == 'FP']),
            len(df[df['confusion_matrix'] == 'FN'])]  # keep this order as in the config.json color_dict
ft.pie_simple(counts, colors, '', output)
