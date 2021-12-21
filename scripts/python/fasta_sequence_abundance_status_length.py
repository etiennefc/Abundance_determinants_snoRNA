#!/usr/bin/python3
import pandas as pd
import itertools
""" Split fastas file of expressed or not expressed C/D box snoRNAs according
    to their length (small or long, i.e < or >= 200 nt). """

df = pd.read_csv(snakemake.input.all_features_labels_df, sep='\t')
df = df[df['sno_type'] == 'C/D']
sno_fasta = snakemake.input.sno_fasta
outputs = [snakemake.output.small_expressed_cd, snakemake.output.long_expressed_cd,
            snakemake.output.small_not_expressed_cd, snakemake.output.long_not_expressed_cd]

# Get gene_id of C/D snoRNAs per length and abundance_status
small_expressed_cd = df[(df['sno_length'] < 200) & (df['abundance_cutoff_2'] == 'expressed')].gene_id_sno.to_list()
long_expressed_cd = df[(df['sno_length'] >= 200) & (df['abundance_cutoff_2'] == 'expressed')].gene_id_sno.to_list()
small_not_expressed_cd = df[(df['sno_length'] < 200) & (df['abundance_cutoff_2'] == 'not_expressed')].gene_id_sno.to_list()
long_not_expressed_cd = df[(df['sno_length'] >= 200) & (df['abundance_cutoff_2'] == 'not_expressed')].gene_id_sno.to_list()

# Create dict from sno_fasta as key: val --> id: sequence
d = {}
with open(sno_fasta, 'r') as f:
    for fasta_id, sequence in itertools.zip_longest(*[f]*2):
            fasta_id = fasta_id.lstrip('>').rstrip('\n')
            sequence = sequence.rstrip('\n')
            d[fasta_id] = sequence

# Get sequence for all snoRNA groups in their respective output file
for i, group in enumerate([small_expressed_cd, long_expressed_cd, small_not_expressed_cd, long_not_expressed_cd]):
    with open(outputs[i], 'w') as f:
        for j, id in enumerate(group):
            seq = d[id]
            f.write(f'>{id}\n')
            f.write(f'{seq}\n')
