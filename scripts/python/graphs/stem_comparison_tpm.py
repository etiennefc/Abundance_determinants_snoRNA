#!/usr/bin/python3
import pandas as pd
import functions as ft
from scipy import stats as st

tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
feature_df = pd.read_csv(snakemake.input.all_features, sep='\t')

# Create average TPM column in tpm_df and merge that column to feature_df
tpm_df['avg_tpm'] = tpm_df.filter(regex='^[A-Z].*_[1-3]$', axis=1).mean(axis=1)
tpm_df = tpm_df[['gene_id', 'avg_tpm']]
feature_df = feature_df.merge(tpm_df, how='left', left_on='gene_id_sno', right_on='gene_id')
feature_df = feature_df.drop(['gene_id'], axis=1)

# Get expressed H/ACA and C/D snoRNAs
haca = feature_df[(feature_df['abundance_cutoff_2'] == "expressed") & (feature_df['sno_type'] == "H/ACA")]
cd = feature_df[(feature_df['abundance_cutoff_2'] == "expressed") & (feature_df['sno_type'] == "C/D")]

# Split according to terminal stem MFE (if < (strong stem) or >= (weak) average MFE)
avg_stem_mfe_haca, avg_stem_mfe_cd = haca['terminal_stem_mfe'].mean(), cd['terminal_stem_mfe'].mean()
haca.loc[haca['terminal_stem_mfe'] < avg_stem_mfe_haca, 'terminal_stem_strength'] = 'Strong'
haca.loc[haca['terminal_stem_mfe'] >= avg_stem_mfe_haca, 'terminal_stem_strength'] = 'Weak'
cd.loc[cd['terminal_stem_mfe'] < avg_stem_mfe_cd, 'terminal_stem_strength'] = 'Strong'
cd.loc[cd['terminal_stem_mfe'] >= avg_stem_mfe_cd, 'terminal_stem_strength'] = 'Weak'
stem_haca, no_stem_haca = haca[haca['terminal_stem_strength'] == 'Strong'], haca[haca['terminal_stem_strength'] == 'Weak']
stem_cd, no_stem_cd = cd[cd['terminal_stem_strength'] == 'Strong'], cd[cd['terminal_stem_strength'] == 'Weak']

# Create violin plots to compare weak and strong terminal stem snoRNA's abundance per snoRNA type
ft.violin(haca, "terminal_stem_strength", "avg_tpm", None, None, "Terminal stem strength",
                "Average abundance across tissues (TPM)", "Abundance of H/ACA snoRNAs with strong and weak terminal stem",
                None, None, snakemake.output.violin_haca)

ft.violin(cd, "terminal_stem_strength", "avg_tpm", None, None, "Terminal stem strength",
                "Average abundance across tissues (TPM)", "Abundance of C/D snoRNAs with strong and weak terminal stem",
                None, None, snakemake.output.violin_cd)

# Compute the significance between groups of H/ACA with a strong or weak terminal stem (same with C/D)
MW_U_stats, p_val = st.mannwhitneyu(stem_haca['avg_tpm'], no_stem_haca['avg_tpm'])
print('H/ACA')
print('Mann-Whitney U statistics:', MW_U_stats, ',  p-value:', p_val)

MW_U_stats, p_val = st.mannwhitneyu(stem_cd['avg_tpm'], no_stem_cd['avg_tpm'])
print('C/D')
print('Mann-Whitney U statistics:', MW_U_stats, ',  p-value:', p_val)
