#!/usr/bin/python3
import pandas as pd

""" Convert T into U in snoRNA sequences and create a fasta file containing all
    snoRNA sequences with their gene_id as names"""

sno_df = pd.read_csv(snakemake.input.snodb, sep='\t')
sno_df = sno_df[['gene_id_sno', 'seq']]

# Replace T by U in snoRNA sequences
sno_df['seq'] = sno_df['seq'].str.replace('T', 'U')

# Create the fasta
dictio = sno_df.set_index('gene_id_sno')['seq'].to_dict()
dictio = {'>'+ k: v for k, v in dictio.items()}  # Add '>' in front of all sno id

with open(snakemake.output.sno_sequences, "a+") as file:  # a+ for append in new file
    for k, v in dictio.items():
        file.write(k+'\n'+v+'\n')
