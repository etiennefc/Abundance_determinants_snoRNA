#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Add an abundance cutoff column for snoRNAs to define them as expressed or not
    expressed in mouse. This column will serve as the label used by the predictor.
    Also find the snoRNA length. """

df = pd.read_csv(snakemake.input.sno_tpm_df, sep='\t')
sno_bed = pd.read_csv(snakemake.input.sno_bed, sep='\t',
                    names=['chr', 'start', 'end', 'gene_id', 'dot',
                            'strand', 'feature', 'dot2', 'attributes'], index_col=False)

# Compute the snoRNA length
sno_bed['sno_length'] = sno_bed['end'].astype(int) - sno_bed['start'].astype(int) + 1
sno_bed = sno_bed[['gene_id', 'sno_length']]
df = df.merge(sno_bed, how='left', on='gene_id')


# If snoRNA is expressed >1 TPM in at least one average sample (average of the triplicates), it is expressed, else not expressed
sno_abundance = df.filter(regex='_[123]$')  # tpm column must end with '_1, _2 or _3'
cols = list(sno_abundance.columns)
triplicates = [cols[n:n+3] for n in range(0, len(cols), 3)]

df['abundance_cutoff'] = ''
for i in range(0, len(df)):
    row = sno_abundance.iloc[i]
    for j, triplicate in enumerate(triplicates):
        if (df.loc[i, triplicate].mean() > 1) == True:
            df.loc[i, 'abundance_cutoff'] = 'expressed'
            break
df['abundance_cutoff'] = df['abundance_cutoff'].replace('', 'not_expressed')

print('Abundance cutoff based on >1 TPM in at least one average condition:')
print(coll.Counter(df['abundance_cutoff']))


df.to_csv(snakemake.output.tpm_label_df, index=False, sep='\t')
