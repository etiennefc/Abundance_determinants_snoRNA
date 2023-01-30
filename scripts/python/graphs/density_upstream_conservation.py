#!/usr/bin/python3
import pandas as pd
import functions as ft
from scipy.stats import mannwhitneyu

sno_cons = pd.read_csv(snakemake.input.sno_cons, sep='\t')
sno_cons = sno_cons.rename(columns={'conservation_score': 'sno_conservation_score'})
upstream_sno_cons = pd.read_csv(snakemake.input.upstream_sno_cons, sep='\t')
upstream_sno_cons = upstream_sno_cons.rename(columns={'conservation_score': 'upstream_conservation_score'})
df = pd.read_csv(snakemake.input.feature_label_df, sep='\t')

# Select intergenic and expressed snoRNAs
#intergenic = df[df['abundance_cutoff_host'] == 'intergenic']
#intergenic = df[(df['abundance_cutoff_host'] == 'intergenic') & (df['abundance_cutoff_2'] == 'not_expressed')]
intergenic = df[(df['abundance_cutoff_host'] == 'intergenic') & (df['abundance_cutoff_2'] == 'expressed')]


# Merge conservation info to intergenic df
cons = sno_cons.merge(upstream_sno_cons, how='left', on='gene_id')
intergenic = intergenic.merge(cons, how='left', left_on='gene_id_sno', right_on='gene_id')

# Separate conserved form recently copied snoRNAs (threshold of 0.5 of sno_conservation_score)
new = intergenic[intergenic['sno_conservation_score'] <= 0.5].upstream_conservation_score
old = intergenic[intergenic['sno_conservation_score'] > 0.5].upstream_conservation_score

# Create density
print(new.median())
print(old.median())

U, pval = mannwhitneyu(list(new), list(old))
print(pval)
nb_old, nb_new = str(len(old)), str(len(new))
colors, groups = ['lightgreen', 'grey'], [f'Conserved snoRNAs\n(conservation score > 0.5) (n={nb_old})', f'Recent snoRNAs\n(conservation score <= 0.5) (n={nb_new})']
ft.density_x([old, new], 'Promoter conservation of expresssed \nintergenic snoRNAs (Conservation score)', 'Density', 'linear', '', colors, groups, snakemake.output.density)
