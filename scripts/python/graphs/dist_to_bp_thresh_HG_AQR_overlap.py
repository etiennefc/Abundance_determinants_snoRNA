#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import collections as coll
from scipy.stats import fisher_exact

cols = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature',
        'dot2', 'gene_info', 'chr_aqr', 'start_aqr', 'end_aqr', 'score_aqr', 'signalValue_aqr', 'strand_aqr', 'pval_aqr']
aqr_df = pd.read_csv(snakemake.input.aqr_overlap_HG, sep='\t', names=cols)
df = pd.read_csv(snakemake.input.df, sep='\t')
exp = df[(df['abundance_cutoff_2'] == 'expressed') & (df['abundance_cutoff_host'] != 'intergenic')]

# Define if sno intron is bound by AQR
aqr_bound_sno_intron = list(pd.unique(aqr_df.gene_id))
exp.loc[exp['gene_id_sno'].isin(aqr_bound_sno_intron), 'AQR_binding'] = 'Intron is bound by AQR'
exp['AQR_binding'] = exp['AQR_binding'].fillna('Intron is not bound by AQR')
cd, haca = exp[exp['sno_type'] == 'C/D'], exp[exp['sno_type'] == 'H/ACA']

# Split expressed snoRNAs based on their distance to bp (close (<=100nt) vs far (>100nt))
cd_close, cd_far = cd[cd['dist_to_bp'] <= 100], cd[cd['dist_to_bp'] > 100]
haca_close, haca_far = haca[haca['dist_to_bp'] <= 100], haca[haca['dist_to_bp'] > 100]



# Generate AQR_binding bar chart comparison between snoRNAs close vs far from bp
cd_close['bp_proximity'] = 'Close to branch\npoint (<= 100 nt)'
cd_far['bp_proximity'] = 'Far from branch\npoint (>100 nt)'
cd_combined = pd.concat([cd_close, cd_far])

counts_per_feature = ft.count_list_x(cd_combined, 'bp_proximity',
                    list(snakemake.params.AQR_binding_colors.keys()),
                    'AQR_binding')
percent = ft.percent_count(counts_per_feature)
ft.stacked_bar(percent, sorted(list(cd_combined['bp_proximity'].unique())),
                list(snakemake.params.AQR_binding_colors.keys()), 'C/D', 'Branch point proximity',
                'Proportion of expressed snoRNAs (%)', snakemake.params.AQR_binding_colors, snakemake.output.bar_HG_cd)


haca_close['bp_proximity'] = 'Close to branch\npoint (<= 100 nt)'
haca_far['bp_proximity'] = 'Far from branch\npoint (>100 nt)'
haca_combined = pd.concat([haca_close, haca_far])
counts_per_feature = ft.count_list_x(haca_combined, 'bp_proximity',
                    list(snakemake.params.AQR_binding_colors.keys()),
                    'AQR_binding')
percent = ft.percent_count(counts_per_feature)
ft.stacked_bar(percent, sorted(list(haca_combined['bp_proximity'].unique())),
                list(snakemake.params.AQR_binding_colors.keys()), 'H/ACA', 'Branch point proximity',
                'Proportion of expressed snoRNAs (%)', snakemake.params.AQR_binding_colors, snakemake.output.bar_HG_haca)



def criteria_count(group1, group2, col, crit):
    "This creates the contingency table needed to perform Fisher's exact test"
    count1a = len(group1[group1[col] == crit])
    count1b = len(group1[group1[col] != crit])
    count2a = len(group2[group2[col] == crit])
    count2b = len(group2[group2[col] != crit])

    dict = {'group1': [count1a, count1b], 'group2': [count2a, count2b]}
    table = pd.DataFrame(data=dict, index=[crit, '!= '+crit])
    print(table)
    return table


for binding_status in list(snakemake.params.AQR_binding_colors.keys()):
    table_cd = criteria_count(cd_close, cd_far, 'AQR_binding', binding_status)
    table_haca = criteria_count(haca_close, haca_far, 'AQR_binding', binding_status)
    oddsratio_cd, p_val_cd = fisher_exact(table_cd)
    oddsratio_haca, p_val_haca = fisher_exact(table_haca)
    print(binding_status)
    print(f"p = {p_val_cd} (Fisher's exact test) C/D")
    print(f"p = {p_val_haca} (Fisher's exact test) H/ACA")



