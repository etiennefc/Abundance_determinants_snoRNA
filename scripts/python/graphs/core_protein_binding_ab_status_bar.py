#!/usr/bin/python3
import functions as ft
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import collections as coll
from scipy.stats import fisher_exact

cols_eclip = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature',
        'dot2', 'gene_info', 'chr_rbp', 'start_rbp', 'end_rbp', 'score_rbp', 'signalValue_rbp', 'strand_rbp', 'pval_rbp']
cols_par_clip = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature',
        'dot2', 'gene_info', 'chr_rbp', 'start_rbp', 'end_rbp', 'strand_rbp', 'score_rbp', 'strand_rbp_duplicated']
dkc1_eclip = pd.read_csv(snakemake.input.dkc1_eclip_overlap, sep='\t', names=cols_eclip)
dkc1_par_clip = pd.read_csv(snakemake.input.dkc1_par_clip_overlap, sep='\t', names=cols_par_clip)
nop58_par_clip = pd.read_csv(snakemake.input.nop58_par_clip_overlap, sep='\t', names=cols_par_clip)
fbl_par_clip = pd.read_csv(snakemake.input.fbl_par_clip_overlap, sep='\t', names=cols_par_clip)
nop56_par_clip = pd.read_csv(snakemake.input.nop56_par_clip_overlap, sep='\t', names=cols_par_clip)
df = pd.read_csv(snakemake.input.df, sep='\t')
cd, haca = df[df['sno_type'] == 'C/D'], df[df['sno_type'] == 'H/ACA']


# Filter PAR-CLIP peaks based on enrichment score
dkc1_par_clip = dkc1_par_clip[dkc1_par_clip['score_rbp'] > 10]
nop58_par_clip = nop58_par_clip[nop58_par_clip['score_rbp'] > 10]
nop56_par_clip = nop56_par_clip[nop56_par_clip['score_rbp'] > 10]
fbl_par_clip = fbl_par_clip[fbl_par_clip['score_rbp'] > 10]


# Define if sno is bound by RBP in given dataset
haca.loc[haca['gene_id_sno'].isin(list(dkc1_eclip.gene_id)), 'DKC1_eCLIP'] = "RBP is bound to snoRNA"
haca.loc[~haca['gene_id_sno'].isin(list(dkc1_eclip.gene_id)), 'DKC1_eCLIP'] = "RBP is not bound to snoRNA"
haca.loc[haca['gene_id_sno'].isin(list(dkc1_par_clip.gene_id)), 'DKC1_PAR_CLIP'] = "RBP is bound to snoRNA"
haca.loc[~haca['gene_id_sno'].isin(list(dkc1_par_clip.gene_id)), 'DKC1_PAR_CLIP'] = "RBP is not bound to snoRNA"

cd.loc[cd['gene_id_sno'].isin(list(nop58_par_clip.gene_id)), 'NOP58_PAR_CLIP'] = "RBP is bound to snoRNA"
cd.loc[~cd['gene_id_sno'].isin(list(nop58_par_clip.gene_id)), 'NOP58_PAR_CLIP'] = "RBP is not bound to snoRNA"
cd.loc[cd['gene_id_sno'].isin(list(fbl_par_clip.gene_id)), 'FBL_PAR_CLIP'] = "RBP is bound to snoRNA"
cd.loc[~cd['gene_id_sno'].isin(list(fbl_par_clip.gene_id)), 'FBL_PAR_CLIP'] = "RBP is not bound to snoRNA"
cd.loc[cd['gene_id_sno'].isin(list(nop56_par_clip.gene_id)), 'NOP56_PAR_CLIP'] = "RBP is bound to snoRNA"
cd.loc[~cd['gene_id_sno'].isin(list(nop56_par_clip.gene_id)), 'NOP56_PAR_CLIP'] = "RBP is not bound to snoRNA"


# Generate RBP_binding bar chart comparison between expressed and not expressed snoRNAs 
def global_stacked_bar(df, hue_col, color_dict, sep_col, title, xlabel, ylabel, output_path):
    counts_per_feature = ft.count_list_x(df, hue_col, list(color_dict.keys()), sep_col)
    percent = ft.percent_count(counts_per_feature)
    ft.stacked_bar(percent, sorted(list(haca[hue_col].unique())),
                    list(color_dict.keys()), title, xlabel,
                    ylabel, color_dict, output_path)

global_stacked_bar(haca, 'abundance_cutoff_2', snakemake.params.RBP_binding_colors, 'DKC1_eCLIP', 'H/ACA with DKC1 eCLIP', 
                    'Expression status of snoRNAs', 'Proportion of expressed snoRNAs (%)', snakemake.output.bar_haca_dkc1_eclip)

global_stacked_bar(haca, 'abundance_cutoff_2', snakemake.params.RBP_binding_colors, 'DKC1_PAR_CLIP', 'H/ACA with DKC1 PAR-CLIP',
                    'Expression status of snoRNAs', 'Proportion of expressed snoRNAs (%)', snakemake.output.bar_haca_dkc1_par_clip)

global_stacked_bar(cd, 'abundance_cutoff_2', snakemake.params.RBP_binding_colors, 'NOP58_PAR_CLIP', 'C/D with NOP58 PAR-CLIP',
                    'Expression status of snoRNAs', 'Proportion of expressed snoRNAs (%)', snakemake.output.bar_cd_nop58_par_clip)

global_stacked_bar(cd, 'abundance_cutoff_2', snakemake.params.RBP_binding_colors, 'FBL_PAR_CLIP', 'C/D with FBL PAR-CLIP',
                    'Expression status of snoRNAs', 'Proportion of expressed snoRNAs (%)', snakemake.output.bar_cd_fbl_par_clip)

global_stacked_bar(cd, 'abundance_cutoff_2', snakemake.params.RBP_binding_colors, 'NOP56_PAR_CLIP', 'C/D with NOP56 PAR-CLIP',
                    'Expression status of snoRNAs', 'Proportion of expressed snoRNAs (%)', snakemake.output.bar_cd_nop56_par_clip)



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


for binding_status in list(snakemake.params.RBP_binding_colors.keys()):
    table_haca_eclip = criteria_count(haca[haca['abundance_cutoff_2'] == 'expressed'], haca[haca['abundance_cutoff_2'] == 'not_expressed'], 'DKC1_eCLIP', binding_status)
    table_haca_par_clip = criteria_count(haca[haca['abundance_cutoff_2'] == 'expressed'], haca[haca['abundance_cutoff_2'] == 'not_expressed'], 'DKC1_PAR_CLIP', binding_status)
    table_cd_nop58 = criteria_count(cd[cd['abundance_cutoff_2'] == 'expressed'], cd[cd['abundance_cutoff_2'] == 'not_expressed'], 'NOP58_PAR_CLIP', binding_status)
    table_cd_fbl = criteria_count(cd[cd['abundance_cutoff_2'] == 'expressed'], cd[cd['abundance_cutoff_2'] == 'not_expressed'], 'FBL_PAR_CLIP', binding_status)
    table_cd_nop56 = criteria_count(cd[cd['abundance_cutoff_2'] == 'expressed'], cd[cd['abundance_cutoff_2'] == 'not_expressed'], 'NOP56_PAR_CLIP', binding_status)

    oddsratio_haca_eclip, p_val_haca_eclip = fisher_exact(table_haca_eclip)
    oddsratio_haca_par_clip, p_val_haca_par_clip = fisher_exact(table_haca_par_clip)
    oddsratio_cd_nop58, p_val_cd_nop58 = fisher_exact(table_cd_nop58)
    oddsratio_cd_fbl, p_val_cd_fbl = fisher_exact(table_cd_fbl)
    oddsratio_cd_nop56, p_val_cd_nop56 = fisher_exact(table_cd_nop56)

    print('\n'+binding_status)
    print(f"p = {p_val_haca_eclip} (Fisher's exact test) H/ACA eCLIP")
    print(f"p = {p_val_haca_par_clip} (Fisher's exact test) H/ACA PAR-CLIP")
    print(f"p = {p_val_cd_nop58} (Fisher's exact test) C/D NOP58 PAR-CLIP")
    print(f"p = {p_val_cd_fbl} (Fisher's exact test) C/D FBL PAR-CLIP")
    print(f"p = {p_val_cd_nop56} (Fisher's exact test) C/D NOP56 PAR-CLIP")


