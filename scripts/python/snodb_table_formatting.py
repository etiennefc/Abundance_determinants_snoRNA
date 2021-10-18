#!/usr/bin/python3
import pandas as pd

""" Drop NA and duplicates from snodb table and format the rest of the columns.
    Use the final host gene list v101 to determine snoRNA host gene biotype,
    function, name, chr, start and end. Add also information about NMD
    susceptibility of host genes and presence of dual-initiation (DI) promoters
    within host genes. Of note, this table contains all snoRNAs and also some
    scaRNAs found on the database snoDB."""

snodb_original = pd.read_csv(snakemake.input.snodb_original_table, sep='\t')
hg_df = pd.read_csv(snakemake.input.host_gene_df)
nmd = pd.read_csv(snakemake.input.nmd)
di_promoter = pd.read_csv(snakemake.input.di_promoter)
hg_functions = pd.read_csv(snakemake.input.host_functions, sep='\t')


# Remove NA lines and duplicates
snodb = snodb_original.dropna(subset=['gene_id_annot2020', 'gene_name_annot2020'])
snodb = snodb.drop_duplicates(subset=['gene_id_annot2020', 'gene_name_annot2020'])


# Drop unwanted columns and rename other columns
snodb = snodb[['gene_id_annot2020', 'gene_name_annot2020', 'box type',
                'target summary', 'seq']]
snodb.columns = ['gene_id_sno', 'gene_name_sno', 'sno_type', 'sno_target_old', 'seq']


# Format the sno_target column
rRNA = ['Others, rRNA', 'rRNA', 'Others, rRNA, snRNA', 'rRNA, snRNA']
snRNA = ['Others, snRNA', 'snRNA']
snodb.loc[snodb['sno_target_old'].isin(rRNA), 'sno_target'] = 'rRNA'
snodb.loc[snodb['sno_target_old'].isin(snRNA), 'sno_target'] = 'snRNA'
snodb.loc[snodb['sno_target_old'] == 'Others', 'sno_target'] = 'Orphan'
snodb['sno_target'] = snodb['sno_target'].fillna('Orphan')


# Add host gene info (id, name, chr (seqname), strand, start, end, biotype)
snodb = snodb.merge(hg_df, how='left', left_on='gene_id_sno', right_on='sno_id')

non_coding = ['lncRNA', 'transcribed_unprocessed_pseudogene', 'unprocessed_pseudogene',
                'transcribed_unitary_pseudogene', 'TEC', 'transcribed_processed_pseudogene',
                'processed_pseudogene']
snodb.loc[snodb['host_biotype'] == 'protein_coding', 'host_biotype2'] = 'protein_coding'
snodb.loc[snodb['host_biotype'].isin(non_coding), 'host_biotype2'] = 'non_coding'
snodb['host_biotype2'] = snodb['host_biotype2'].fillna('intergenic')


# Add NMD and DI promoter info for intronic snoRNAs
nmd_substrates = list(pd.unique(nmd['gs']))
di_promoter_host = list(pd.unique(di_promoter['gene_id']))
snodb.loc[snodb['host_name'].isin(nmd_substrates), 'NMD_susceptibility'] = True
snodb.loc[(~snodb['host_name'].isin(nmd_substrates)) &
        (~snodb['host_name'].isnull()), 'NMD_susceptibility'] = False
snodb.loc[snodb['host_name'].isnull(), 'NMD_susceptibility'] = "intergenic"

snodb.loc[snodb['host_id'].isin(di_promoter_host), 'di_promoter'] = "dual_initiation"
snodb.loc[(~snodb['host_id'].isin(di_promoter_host)) &
        (~snodb['host_id'].isnull()), 'di_promoter'] = "simple_initiation"
snodb.loc[snodb['host_id'].isnull(), 'di_promoter'] = "intergenic"


# Add host_function column (combined for protein-coding and non-coding HGs)
snodb = snodb.merge(hg_functions, how='left', left_on='host_id', right_on='host_id')
snodb['host_function'] = snodb['host_function'].fillna('intergenic')

snodb.drop(columns=['sno_target_old', 'sno_id', 'host_biotype'], inplace=True)

snodb.to_csv(snakemake.output.snodb_formatted, index=False, sep='\t')
