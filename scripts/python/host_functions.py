#!/usr/bin/python3
import pandas as pd

""" Generate a host_function column for all snoRNAs. Protein-coding HGs are
    separated into 5 categories (ribosomal protein; ribosome biogenesis &
    translation; RNA binding, processing & splicing; Other; poorly
    characterized). Non-coding HGs are separated in 2 categories
    (functional_ncRNA; non_functional_ncRNA). For intergenic snoRNAs, they will
    have 'intergenic' as host function"""

host_gene_df = pd.read_csv(snakemake.input.host_genes)
host_gene_df = host_gene_df[['host_id', 'host_biotype']]
non_coding = pd.read_csv(snakemake.input.nc_functions, sep='\t', header=0, names=['nc_functions'])
protein_coding = pd.read_csv(snakemake.input.pc_functions, sep='\t')
protein_coding = protein_coding.drop_duplicates(subset=['host_id'])
protein_coding = protein_coding[['host_id', 'pc_host_function']]

# Merge non-coding functions to host gene df
host_gene_df = host_gene_df.merge(non_coding, how='left', left_on='host_id',
                                    right_on='nc_functions')

# Merge protein-coding functions to host gene df
host_gene_df = host_gene_df.merge(protein_coding, how='left', left_on='host_id',
                                    right_on='host_id')

# Create a host_function column
host_gene_df.loc[(host_gene_df['host_biotype'] != 'protein_coding') &
                (host_gene_df['nc_functions'].notnull()), 'host_function'] = 'functional_ncRNA'
host_gene_df.loc[(host_gene_df['host_biotype'] != 'protein_coding') &
                (host_gene_df['nc_functions'].isnull()), 'host_function'] = 'non_functional_ncRNA'

host_gene_df.loc[host_gene_df['pc_host_function'] == 'Ribosomal protein',
                'host_function'] = 'Ribosomal protein'
host_gene_df.loc[host_gene_df['pc_host_function'] == 'Ribosome biogenesis & translation',
                'host_function'] = 'Ribosome biogenesis & translation'
host_gene_df.loc[host_gene_df['pc_host_function'] == 'RNA binding, processing, splicing',
                'host_function'] = 'RNA binding, processing, splicing'
host_gene_df.loc[host_gene_df['pc_host_function'] == 'Other', 'host_function'] = 'Other'
host_gene_df.loc[host_gene_df['pc_host_function'] == 'Poorly characterized',
                'host_function'] = 'Poorly characterized'

host_gene_df = host_gene_df[['host_id', 'host_function']].drop_duplicates(subset=['host_id'])

host_gene_df.to_csv(snakemake.output.hg_function_df, index=False, sep='\t')
