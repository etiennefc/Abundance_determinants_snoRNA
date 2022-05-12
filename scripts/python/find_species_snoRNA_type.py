#!/usr/bin/python3
import pandas as pd

""" Find snoRNA type (C/D vs H/ACA) of species snoRNAs """

# Loading file and dfs
rna_central_df = pd.read_csv(snakemake.input.rna_central_df, sep='\t')
id_conversion_df = pd.read_csv(snakemake.input.id_conversion_df, sep='\t',
                        names=['rna_central_id', 'source', 'ensembl_transcript_id',
                                'taxid', 'biotype', 'ensembl_gene_id'])
sno_bed = pd.read_csv(snakemake.input.sno_bed, sep='\t',
                        names=['chr', 'start', 'end', 'gene_id', 'dot', 'strand',
                                'source', 'feature', 'dot2', 'attributes'])

attributes = [att.split(';') for att in list(sno_bed.attributes)]
gene_id_name = {}
for sno_attributes in attributes:
    temp_name = None
    for specific_attribute in sno_attributes:
        if 'gene_id' in specific_attribute:
            id = specific_attribute.split(' "')[-1].strip('"')
        elif 'gene_name' in specific_attribute:
            name = specific_attribute.split(' "')[-1].strip('"')
            temp_name = name
    if all('gene_name' not in attri for attri in sno_attributes):
        temp_name = id
    gene_id_name[id] = temp_name


# Get snoRNA ids and name and create df
sno_df = pd.DataFrame(gene_id_name.items(), columns=['gene_id', 'gene_name'])

# Remove ".number" a the end of ensembl gene id
if '.' in list(id_conversion_df['ensembl_gene_id'])[0]:
    id_conversion_df[['ensembl_gene_id', 'suffix']] = id_conversion_df['ensembl_gene_id'].str.split('.', expand=True)

# Dict of corresponding ensembl and rnacentral ids
id_mapping = dict(zip(id_conversion_df.ensembl_gene_id, id_conversion_df.rna_central_id))

# Find if sno_type is known for species snoRNAs in RNAcentral description
rna_central_df.loc[rna_central_df['description'].str.contains('C/D|SNORD|Snord|U3'), 'snoRNA_type'] = 'C/D'
rna_central_df.loc[rna_central_df['description'].str.contains('H/ACA|SNORA|Snora|ACA'), 'snoRNA_type'] = 'H/ACA'
sno_type_dict = dict(zip(rna_central_df.upi, rna_central_df.snoRNA_type))  # where upi is the column of RNA central id

# Create rna_central_id and snoRNA_type columns
sno_df['rna_central_id'] = sno_df['gene_id'].map(id_mapping)
sno_df['snoRNA_type'] = sno_df['rna_central_id'].map(sno_type_dict)

# For snoRNAs without snoRNA type from RNAcentral description, find if their gene_name contains this information directly
sno_df.loc[(sno_df['snoRNA_type'].isna()) & (sno_df['gene_name'].str.contains('Snord|SNORD|U8|U3')), 'snoRNA_type'] = 'C/D'
sno_df.loc[(sno_df['snoRNA_type'].isna()) & (sno_df['gene_name'].str.contains('Snora|SNORA')), 'snoRNA_type'] = 'H/ACA'


# Drop all snoRNAs where we couldn't find the snoRNA type
sno_df = sno_df.dropna(subset=['snoRNA_type'])


sno_df.to_csv(snakemake.output.snoRNA_type_df, sep='\t', index=False)
