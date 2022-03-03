#!/usr/bin/python3
import pandas as pd
import functions as ft

df = pd.read_csv(snakemake.input.df, sep='\t')
feature = snakemake.wildcards.intronic_features
limits = [[0, 25, '[0; 25['], [25, 50, '[25; 50['], [50, 100, '[50; 100['],
            [100, 250, '[100; 250['], [250, 500, '[250; 500['],
            [500, 1000, '[500; 1000['], [1000, 5000, '[1000; 5000['],
            [5000, 10000, '[5000; 10000['], [10000, 100000, '[10000; 100000['],
            [100000, 1000000000, '>100000']]
labels = []
# Create ranges of data for intronic features (with large scale of data)
for i, limit in enumerate(limits):
    lower, upper, label = limit[0], limit[1], limit[2]
    df.loc[(df[feature] >= lower) & (df[feature] < upper), feature+'_range'] = label
    labels.append(label)

# Drop intergenic snoRNAs
df = df.dropna(subset=[feature])
cd = df[df['sno_type'] == 'C/D']
haca = df[df['sno_type'] == 'H/ACA']


# Generate a grouped bar chart of intronic features with a hue of abundance_cutoff_2
# for either C/D and H/ACA snoRNAs separately

#For C/D
counts_cd = ft.count_list_x_unsorted(cd, list(snakemake.params.hue_color.keys()),
                'abundance_cutoff_2', labels, feature+'_range')

ft.bar_from_lst_of_lst(counts_cd, [0.15, 0.45], list(snakemake.params.hue_color.values()), 0.3, labels,
                        feature+'_range', "Number of snoRNAs",
                        list(snakemake.params.hue_color.keys()), snakemake.output.bar_cd)


#For H/ACA
counts_haca = ft.count_list_x_unsorted(haca, list(snakemake.params.hue_color.keys()),
                'abundance_cutoff_2', labels, feature+'_range')

ft.bar_from_lst_of_lst(counts_haca, [0.15, 0.45], list(snakemake.params.hue_color.values()), 0.3, labels,
                        feature+'_range', "Number of snoRNAs",
                        list(snakemake.params.hue_color.keys()), snakemake.output.bar_haca)
