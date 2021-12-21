#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Add an abundance cutoff column for all RNAs to define them as expressed or not
    expressed. """

df = pd.read_csv(snakemake.input.tpm_df)

# If oRNA is expressed >1 TPM in at least one average tissue (average of the triplicates), it is expressed, else not expressed
rna_abundance = df.iloc[:, 4:25]
cols = list(rna_abundance.columns)
triplicates = [cols[n:n+3] for n in range(0, len(cols), 3)]

df['abundance_cutoff_2'] = ''
for i in range(0, len(df)):
    row = rna_abundance.iloc[i]
    for j, triplicate in enumerate(triplicates):
        if (df.loc[i, triplicate].mean() > 1) == True:
            df.loc[i, 'abundance_cutoff_2'] = 'expressed'
            break
df['abundance_cutoff_2'] = df['abundance_cutoff_2'].replace('', 'not_expressed')

for i, bio in enumerate(list(pd.unique(df['gene_biotype']))):
    print(bio)
    d = df[df['gene_biotype'] == bio]
    print(coll.Counter(d['abundance_cutoff_2']))  


df.to_csv(snakemake.output.abundance_cutoff_df, index=False, sep='\t')
