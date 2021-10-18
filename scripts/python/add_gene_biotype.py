#!/usr/bin/python3
import pandas as pd


"""Add a gene biotype and a simplified gene biotype column in an existing dataframe using a reference table"""

ref_table = pd.read_csv(snakemake.input.ref_table, sep='\t')
df = pd.read_csv(snakemake.input.tpm_df)

# Create a dictionary of gene_id/gene_biotype (key/value) from the ref table and use it to create the gene_biotype col
ref_dict = ref_table.set_index('gene_id').to_dict()['gene_biotype']
df.insert(loc=2, column='gene_biotype', value=df['gene_id'].map(ref_dict))  # insert new col at position 2

# Create a dictionary of simplified gene biotypes and use it to create the gene_biotype2 column
simplified_dict = {}
all_biotypes = list(set(ref_dict.values()))
pc_list = ['IG_C_gene', 'IG_D_gene', 'TR_V_gene', 'IG_V_gene', 'IG_J_gene', 'TR_J_gene', 'TR_C_gene', 'TR_D_gene', 'protein_coding', ]
pseudogene_list = ['translated_unprocessed_pseudogene', 'translated_processed_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene',
                   'processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'transcribed_unitary_pseudogene', 'transcribed_processed_pseudogene',
                   'IG_V_pseudogene', 'pseudogene', 'TR_V_pseudogene', 'TR_J_pseudogene', 'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene',
                   'polymorphic_pseudogene', 'rRNA_pseudogene']
other_list = ['rRNA', 'Mt_rRNA', 'ribozyme', 'scRNA', 'vault_RNA', 'sRNA', 'ETS-RNA', 'ITS-RNA', 'TEC']
tRNA_list = ['tRNA', 'Mt_tRNA', 'pre-tRNA', 'tRNA_fragment']
same_biotype = ['lncRNA', 'misc_RNA', 'intronic_cluster', 'intergenic_cluster', 'snoRNA', 'snRNA', 'miRNA', 'scaRNA']


for gene in pc_list:
    simplified_dict[gene] = 'protein_coding'
for gene in pseudogene_list:
    simplified_dict[gene] = 'pseudogene'
for gene in other_list:
    simplified_dict[gene] = 'other'
for gene in tRNA_list:
    simplified_dict[gene] = 'tRNA'
for gene in same_biotype:
    simplified_dict[gene] = gene

df.insert(loc=3, column='gene_biotype2', value=df['gene_biotype'].map(simplified_dict))  # insert new col at position 3

df.to_csv(snakemake.output.tpm_biotype, index=False)
