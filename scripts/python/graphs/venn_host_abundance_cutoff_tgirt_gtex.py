#!/usr/bin/python3
import pandas as pd
import functions_venn as ft
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Create Venn diagram and compute Simple Matching Coefficient (SMC) (how similar two vectors are (intersection of both vectors / union of both vectors))
# SMC of 0 being totally dissimilar and of 1 being totally similar
gtex_df = pd.read_csv(snakemake.input.gtex_df, sep='\t')
tgirt_df = pd.read_csv(snakemake.input.tgirt_df, sep='\t')

# Drop intergenic snoRNAs
gtex_df = gtex_df[gtex_df['intergenic'] == 0]
tgirt_df = tgirt_df[tgirt_df['intergenic'] == 0]


# Get the number of snoRNA host genes expressed in both or either TGIRT and GTEx datasets
# Get also the number of snoRNA host gene not expressed in both datasets

tgirt_only, gtex_only, both_expressed, both_not_expressed = [], [], [], []
for i, tgirt_val in enumerate(tgirt_df.host_expressed):
    gtex_val = list(gtex_df.host_expressed)[i]
    if (tgirt_val == 1) & (gtex_val == 1):
        both_expressed.append(tgirt_val)
    elif (tgirt_val == 1) & (gtex_val == 0):
        tgirt_only.append(tgirt_val)
    elif (tgirt_val == 0) & (gtex_val == 1):
        gtex_only.append(tgirt_val)
    elif (tgirt_val == 0) & (gtex_val == 0):
        both_not_expressed.append(tgirt_val)


# Compute SMC
smc = (len(both_expressed) + len(both_not_expressed)) / (len(both_expressed) + len(both_not_expressed) + len(tgirt_only) + len(gtex_only))
smc = str(smc)

# Create host_expressed Venn diagram
ft.venn_2([len(tgirt_only), len(gtex_only), len(both_expressed)],
            ['lightblue', 'red'], ['TGIRT', 'GTEx'],
            f'Number of snoRNA host genes expressed\nin TGIRT and GTEx datasets (SMC={smc})',
            snakemake.output.venn_host_expressed)

# Clear the axis between the saving of the second Venn diagram
plt.clf()

# Create host_not_expressed Venn diagram
# (tgirt_only is the opposite of gtex_only and vice-versa (this is why we switch the order compared to the previous venn diagram))
# this means that the number of HG only expressed in TGIRT is equal to the number of HG not expressed only in GTEx and vice-versa
ft.venn_2([len(gtex_only), len(tgirt_only), len(both_not_expressed)],
        ['lightblue', 'red'], ['TGIRT', 'GTEx'],
        f'Number of snoRNA host genes not expressed\nin TGIRT and GTEx datasets (SMC={smc})',
        snakemake.output.venn_host_not_expressed)
