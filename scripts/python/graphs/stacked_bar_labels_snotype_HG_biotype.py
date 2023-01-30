#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import collections as coll

df = pd.read_csv(snakemake.input.df, sep='\t')
df = df[['abundance_cutoff_2', 'sno_type', 'host_biotype2']]
exp, not_exp = df[df['abundance_cutoff_2'] == 'expressed'], df[df['abundance_cutoff_2'] == 'not_expressed']
df_len = len(df)

snotype_dict_exp, snotype_dict_not_exp = {k:v for k,v in sorted(dict(coll.Counter(exp.sno_type)).items())}, {k:v for k,v in sorted(dict(coll.Counter(not_exp.sno_type)).items())}
host_biotype_dict_exp, host_biotype_dict_not_exp = {k:v for k,v in sorted(dict(coll.Counter(exp.host_biotype2)).items())}, {k:v for k,v in sorted(dict(coll.Counter(not_exp.host_biotype2)).items())}

# Get values (and as %) of each hue per expression status
sno_type_exp_nb, sno_type_not_exp_nb = list(snotype_dict_exp.values()), list(snotype_dict_not_exp.values())
host_biotype_exp_nb, host_biotype_not_exp_nb = list(host_biotype_dict_exp.values()), list(host_biotype_dict_not_exp.values())

def get_percent_df(l1, l2, index, cols):
    percent1 = [i * 100 / sum(l1) for i in l1]
    percent2 = [i * 100 / sum(l2) for i in l2]
    df = pd.DataFrame([percent1, percent2], index = index, columns = cols)
    print(df)
    return df
exp_len, not_exp_len = str(len(exp)), str(len(not_exp))
ind = [f'Expressed\n({exp_len})', f'Not expressed\n({not_exp_len})']
d1 = get_percent_df(sno_type_exp_nb, sno_type_not_exp_nb, ind, list(snotype_dict_exp.keys()))
d2 = get_percent_df(host_biotype_exp_nb, host_biotype_not_exp_nb, ind, list(host_biotype_dict_exp.keys()))
df_l = [d2, d1]

# Create a grouped stacked bar chart showing the % of either expressed or not experssed snoRNAs (separate bar) 
# The hue of each bar is either the snoRNA type (C/D or H/ACA) or the host gene biotype
colors = [list(snakemake.params.host_biotype_colors.values()), list(snakemake.params.sno_type_colors.values())]
ft.plot_clustered_stacked2(df_l, colors, 'Expression status of snoRNAs', 'Proportion of snoRNAs (%)', ind, snakemake.output.bar)

