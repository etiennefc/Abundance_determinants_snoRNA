#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import collections as coll
from scipy.stats import mannwhitneyu, fisher_exact

colors = list(snakemake.params.dist_to_bp_group_colors.values())
df = pd.read_csv(snakemake.input.df, sep='\t')
exp = df[(df['abundance_cutoff_2'] == 'expressed') & (df['abundance_cutoff_host'] != 'intergenic')]
cd, haca = exp[exp['sno_type'] == 'C/D'], exp[exp['sno_type'] == 'H/ACA']

# Split expressed snoRNAs based on their distance to bp (close (<=100nt) vs far (>100nt))
cd_close, cd_far = cd[cd['dist_to_bp'] <= 100], cd[cd['dist_to_bp'] > 100]
haca_close, haca_far = haca[haca['dist_to_bp'] <= 100], haca[haca['dist_to_bp'] > 100]

# Generate terminal stem mfe density comparison between snoRNAs close vs far from bp
ft.density_x([cd_close['terminal_stem_mfe'], cd_far['terminal_stem_mfe']], 'Terminal stem stability (kcal/mol)',
            'Density', 'linear', 'Expressed C/D', colors, ['Close to branch point (<= 100 nt)', 'Far from branch point (>100 nt)'], snakemake.output.density_terminal_stem_cd)

ft.density_x([haca_close['terminal_stem_mfe'], haca_far['terminal_stem_mfe']], 'Terminal stem stability (kcal/mol)',
            'Density', 'linear', 'Expressed H/ACA', colors, ['Close to branch point (<= 100 nt)', 'Far from branch point (>100 nt)'], snakemake.output.density_terminal_stem_haca)

U, p = mannwhitneyu(cd_close['terminal_stem_mfe'], cd_far['terminal_stem_mfe'])
print(f'p = {p} (M-W U test, C/D close vs C/D far, terminal stem mfe)')

U, p = mannwhitneyu(haca_close['terminal_stem_mfe'], haca_far['terminal_stem_mfe'])
print(f'p = {p} (M-W U test, H/ACA close vs H/ACA far, terminal stem mfe)')

# Generate box score density comparison between snoRNAs close vs far from bp
ft.density_x([cd_close['combined_box_hamming'], cd_far['combined_box_hamming']], 'Box score',
            'Density', 'linear', 'Expressed C/D', colors, ['Close to branch point (<= 100 nt)', 'Far from branch point (>100 nt)'], snakemake.output.density_box_score_cd)

ft.density_x([haca_close['combined_box_hamming'], haca_far['combined_box_hamming']], 'Box score',
            'Density', 'linear', 'Expressed H/ACA', colors, ['Close to branch point (<= 100 nt)', 'Far from branch point (>100 nt)'], snakemake.output.density_box_score_haca)

U, p = mannwhitneyu(cd_close['combined_box_hamming'], cd_far['combined_box_hamming'])
print(f'p = {p} (M-W U test, C/D close vs C/D far, Box score)')

U, p = mannwhitneyu(haca_close['combined_box_hamming'], haca_far['combined_box_hamming'])
print(f'p = {p} (M-W U test, H/ACA close vs H/ACA far, Box score)')


# Generate sno mfe density comparison between snoRNAs close vs far from bp
ft.density_x([cd_close['sno_mfe'], cd_far['sno_mfe']], 'Structure stability (kcal/mol)',
            'Density', 'linear', 'Expressed C/D', colors, ['Close to branch point (<= 100 nt)', 'Far from branch point (>100 nt)'], snakemake.output.density_sno_mfe_cd)

ft.density_x([haca_close['sno_mfe'], haca_far['sno_mfe']], 'Structure stability (kcal/mol)',
            'Density', 'linear', 'Expressed H/ACA', colors, ['Close to branch point (<= 100 nt)', 'Far from branch point (>100 nt)'], snakemake.output.density_sno_mfe_haca)

U, p = mannwhitneyu(cd_close['sno_mfe'], cd_far['sno_mfe'])
print(f'p = {p} (M-W U test, C/D close vs C/D far, Sno MFE)')

U, p = mannwhitneyu(haca_close['sno_mfe'], haca_far['sno_mfe'])
print(f'p = {p} (M-W U test, H/ACA close vs H/ACA far, Sno MFE)')


# Generate target type bar chart comparison between snoRNAs close vs far from bp
cd_close['bp_proximity'] = 'Close to branch\npoint (<= 100 nt)'
cd_far['bp_proximity'] = 'Far from branch\npoint (>100 nt)'
cd_combined = pd.concat([cd_close, cd_far])
counts_per_feature = ft.count_list_x(cd_combined, 'bp_proximity',
                    list(snakemake.params.sno_target_colors.keys()),
                    'sno_target')
percent = ft.percent_count(counts_per_feature)
ft.stacked_bar(percent, sorted(list(cd_combined['bp_proximity'].unique())),
                list(snakemake.params.sno_target_colors.keys()), 'C/D', 'Branch point proximity',
                'Proportion of expressed snoRNAs (%)', snakemake.params.sno_target_colors, snakemake.output.bar_target_cd)


haca_close['bp_proximity'] = 'Close to branch\npoint (<= 100 nt)'
haca_far['bp_proximity'] = 'Far from branch\npoint (>100 nt)'
haca_combined = pd.concat([haca_close, haca_far])
counts_per_feature = ft.count_list_x(haca_combined, 'bp_proximity',
                    list(snakemake.params.sno_target_colors.keys()),
                    'sno_target')
percent = ft.percent_count(counts_per_feature)
print(percent)

ft.stacked_bar(percent, sorted(list(haca_combined['bp_proximity'].unique())),
                list(snakemake.params.sno_target_colors.keys()), 'H/ACA', 'Branch point proximity',
                'Proportion of expressed snoRNAs (%)', snakemake.params.sno_target_colors, snakemake.output.bar_target_haca)

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


for target in list(snakemake.params.sno_target_colors.keys()):
    table_cd = criteria_count(cd_close, cd_far, 'sno_target', target)
    table_haca = criteria_count(haca_close, haca_far, 'sno_target', target)
    oddsratio_cd, p_val_cd = fisher_exact(table_cd)
    oddsratio_haca, p_val_haca = fisher_exact(table_haca)
    print(target)
    print(f"p = {p_val_cd} (Fisher's exact test) C/D")
    print(f"p = {p_val_haca} (Fisher's exact test) H/ACA")



