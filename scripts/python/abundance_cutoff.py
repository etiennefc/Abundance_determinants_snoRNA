#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Add an abundance cutoff column for snoRNAs to define them as expressed or not
    expressed. This column will serve as the label used by the predictor. Also
    define an abundance cutoff for host genes that will be used as a feature. """

df = pd.read_csv(snakemake.input.tpm_df)

# If snoRNA is expressed >1 TPM in at least one tissue sample, it is expressed, else not expressed
df.loc[(df.iloc[:, 4:25] > 1).any(axis=1), 'abundance_cutoff'] = 'expressed'
df.loc[(df.iloc[:, 4:25] <= 1).all(axis=1), 'abundance_cutoff'] = 'not_expressed'
print('Abundance cutoff based on >1 TPM in at least one sample:')
print(coll.Counter(df['abundance_cutoff']))  # 967 not_expressed; 574 expressed


# If snoRNA is expressed >1 TPM in at least one average tissue (average of the triplicates), it is expressed, else not expressed
sno_abundance = df.iloc[:, 4:25]
cols = list(sno_abundance.columns)
triplicates = [cols[n:n+3] for n in range(0, len(cols), 3)]

df['abundance_cutoff_2'] = ''
for i in range(0, len(df)):
    row = sno_abundance.iloc[i]
    for j, triplicate in enumerate(triplicates):
        if (df.loc[i, triplicate].mean() > 1) == True:
            df.loc[i, 'abundance_cutoff_2'] = 'expressed'
            break
df['abundance_cutoff_2'] = df['abundance_cutoff_2'].replace('', 'not_expressed')

print('Abundance cutoff based on >1 TPM in at least one average tissue:')
print(coll.Counter(df['abundance_cutoff_2']))  # 1056 not expressed; 485 expressed

# If host gene is expressed >1 TPM in at least one average tissue (average of the triplicates), it is expressed, else not expressed (or no host gene at all)
hg_abundance = df.iloc[:, 33:54]
cols_hg = list(hg_abundance.columns)
triplicates_hg = [cols_hg[n:n+3] for n in range(0, len(cols_hg), 3)]

df['abundance_cutoff_host'] = ''
for i in range(0, len(df)):
    row = hg_abundance.iloc[i]
    for j, triplicate in enumerate(triplicates_hg):
        if (df.loc[i, triplicate].mean() > 1) == True:
            df.loc[i, 'abundance_cutoff_host'] = 'host_expressed'
            break
df['abundance_cutoff_host'] = df['abundance_cutoff_host'].replace('', 'host_not_expressed')
df.loc[df['host_id'].isnull(), 'abundance_cutoff_host'] = 'intergenic'


df.to_csv(snakemake.output.abundance_cutoff_df, index=False, sep='\t')
