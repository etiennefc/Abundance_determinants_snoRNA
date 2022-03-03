#!/usr/bin/python3
import pandas as pd
import functions as ft
from scipy.stats import fisher_exact

FP_path = [path for path in snakemake.input.confusion_value_df if 'FP' in path][0]
TN_path = [path for path in snakemake.input.confusion_value_df if 'TN' in path][0]
FP, TN = pd.read_csv(FP_path, sep='\t'), pd.read_csv(TN_path, sep='\t')
multi_HG_df = pd.read_csv(snakemake.input.multi_HG_df, sep='\t')
multi_HG_sno = list(multi_HG_df.gene_id_sno)
color_dict = snakemake.params.color_dict

# Create tables required for the Fisher's exact test contingency table
FP.loc[FP['gene_id_sno'].isin(multi_HG_sno), 'multi_HG_different_labels'] = 'yes'
FP['multi_HG_different_labels'] = FP['multi_HG_different_labels'].fillna('no')

TN.loc[TN['gene_id_sno'].isin(multi_HG_sno), 'multi_HG_different_labels'] = 'yes'
TN['multi_HG_different_labels'] = TN['multi_HG_different_labels'].fillna('no')

# Create contingency table
table = ft.fisher_contingency(FP, TN, 'multi_HG_different_labels', 'yes')
oddsratio, p_val = fisher_exact(table)

print(f'p-val is {p_val}')

# Create bar chart
counts_per_conf_val = [list(table['group1']), list(table['group2'])]
percent = ft.percent_count(counts_per_conf_val)

ft.stacked_bar(percent, ['FP', 'TN'], color_dict.keys(),
                '', 'Confusion value',
                'Proportion of snoRNAs within \n multi-intronic HG (%)', color_dict.values(),
            snakemake.output.bar)
