#!/usr/bin/python3
import pandas as pd
import itertools
""" Extract specific snoRNA sequences (per sno_type and abundance_status) from
    the fasta file of all snoRNA sequences."""

df = pd.read_csv(snakemake.input.all_features_labels_df, sep='\t')
sno_fasta = snakemake.input.sno_fasta
outputs = [snakemake.output.expressed_cd, snakemake.output.expressed_haca,
            snakemake.output.not_expressed_cd, snakemake.output.not_expressed_haca]

# Get gene_id of snoRNAs per sno_type and abundance_status
expressed_cd = df[(df['sno_type'] == 'C/D') & (df['abundance_cutoff_2'] == 'expressed')].gene_id_sno.to_list()
expressed_haca = df[(df['sno_type'] == 'H/ACA') & (df['abundance_cutoff_2'] == 'expressed')].gene_id_sno.to_list()
not_expressed_cd = df[(df['sno_type'] == 'C/D') & (df['abundance_cutoff_2'] == 'not_expressed')].gene_id_sno.to_list()
not_expressed_haca = df[(df['sno_type'] == 'H/ACA') & (df['abundance_cutoff_2'] == 'not_expressed')].gene_id_sno.to_list()

# Create dict from sno_fasta as key: val --> id: sequence
d = {}
with open(sno_fasta, 'r') as f:
    for fasta_id, sequence in itertools.zip_longest(*[f]*2):
            fasta_id = fasta_id.lstrip('>').rstrip('\n')
            sequence = sequence.rstrip('\n')
            d[fasta_id] = sequence

# Get sequence for all snoRNA groups in their respective output file
for i, group in enumerate([expressed_cd, expressed_haca, not_expressed_cd, not_expressed_haca]):
    with open(outputs[i], 'w') as f:
        for j, id in enumerate(group):
            seq = d[id]
            f.write(f'>{id}\n')
            f.write(f'{seq}\n')
