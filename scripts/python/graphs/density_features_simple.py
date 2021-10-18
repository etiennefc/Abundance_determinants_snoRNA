#!/usr/bin/python3
import pandas as pd
import functions as ft
import numpy as np

df = pd.read_csv(snakemake.input.df, sep='\t')

# Generate a simple density plot of numerical_features without a hue
ft.density(df[snakemake.wildcards.numerical_features],
        snakemake.wildcards.numerical_features, 'Density', '',
        snakemake.output.density_features_simple, color=snakemake.params.simple_color)
