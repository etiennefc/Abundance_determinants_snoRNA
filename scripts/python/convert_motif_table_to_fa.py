#!/usr/bin/python3
import pandas as pd

c_d_box_expressed = pd.read_csv(snakemake.input.c_d_box_location_expressed, sep='\t')
c_d_box_not_expressed = pd.read_csv(snakemake.input.c_d_box_location_not_expressed, sep='\t')
h_aca_box_expressed = pd.read_csv(snakemake.input.h_aca_box_location_expressed, sep='\t')
h_aca_box_not_expressed = pd.read_csv(snakemake.input.h_aca_box_location_not_expressed, sep='\t')

def table_to_fasta(df, col1, col2, output_path):
    """ Use 2 columns in table to create fasta (col1 being used for ID lines and
        col2 being sequence lines)."""
    d = dict(zip(df[col1], df[col2]))
    with open(output_path, 'w') as f:
        for sno_id, sequence in d.items():
            f.write(f'>{sno_id}\n')
            f.write(f'{sequence}\n')



table_to_fasta(c_d_box_expressed, 'gene_id', 'C_sequence', snakemake.output.c_expressed)
table_to_fasta(c_d_box_expressed, 'gene_id', 'D_sequence', snakemake.output.d_expressed)
table_to_fasta(c_d_box_not_expressed, 'gene_id', 'C_sequence', snakemake.output.c_not_expressed)
table_to_fasta(c_d_box_not_expressed, 'gene_id', 'D_sequence', snakemake.output.d_not_expressed)
table_to_fasta(c_d_box_expressed, 'gene_id', 'C_prime_sequence', snakemake.output.c_prime_expressed)
table_to_fasta(c_d_box_expressed, 'gene_id', 'D_prime_sequence', snakemake.output.d_prime_expressed)
table_to_fasta(c_d_box_not_expressed, 'gene_id', 'C_prime_sequence', snakemake.output.c_prime_not_expressed)
table_to_fasta(c_d_box_not_expressed, 'gene_id', 'D_prime_sequence', snakemake.output.d_prime_not_expressed)
table_to_fasta(h_aca_box_expressed, 'gene_id', 'H_sequence', snakemake.output.h_expressed)
table_to_fasta(h_aca_box_expressed, 'gene_id', 'ACA_sequence', snakemake.output.aca_expressed)
table_to_fasta(h_aca_box_not_expressed, 'gene_id', 'H_sequence', snakemake.output.h_not_expressed)
table_to_fasta(h_aca_box_not_expressed, 'gene_id', 'ACA_sequence', snakemake.output.aca_not_expressed)
