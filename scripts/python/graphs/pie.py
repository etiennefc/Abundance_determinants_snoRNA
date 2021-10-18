#!/usr/bin/python3
import pandas as pd
import functions as ft
import numpy as np

df = pd.read_csv(snakemake.input.df, sep='\t')
colors = snakemake.params.colors
output = snakemake.output.pie
# Generate a pie chart of the number of snoRNA per abundance status (expressed
# vs not_expressed)
ab_status = [len(df[df['abundance_cutoff_2'] == 'expressed']),
            len(df[df['abundance_cutoff_2'] == 'not_expressed'])]
ft.pie_simple(ab_status, colors, '', output)
