#!/usr/bin/python3
import pandas as pd
import collections as coll
import re
import regex
from math import isnan

""" Find snoRNA type (C/D vs H/ACA) of mouse snoRNAs """

def cut_sequence(seq):
    """ Get the 20 first and 20 last nt of a given sequence."""
    first, last = seq[:20], seq[-20:]
    length = len(seq)
    return first, last, length


def find_c_box(seq):
    """ Find exact C box (RUGAUGA, where R is A or G), if not present, find C
        box with 1,2 or max 3 substitutions. Return also the start and end
        position of that box as 1-based values. If no C box is found, return a
        'NNNNNNN' empty C box and 0 as start and end of C box."""
    first_20, last_20, length_seq = cut_sequence(seq)
    len_c_box = 7
    # First, find exact C box (RUGAUGA) within 20 nt of the snoRNA 5' end
    if re.search('(A|G)UGAUGA', first_20) is not None:  # find exact C box
        i = 1
        for possible_c in re.finditer('(A|G)UGAUGA', first_20):
            if i <= 1:  # select first matched group only (closest RUGAUGA to 5' end of snoRNA)
                c_motif = possible_c.group(0)
                c_start = possible_c.start() + 1
                c_end = possible_c.end()
                i += 1
                return c_motif, c_start, c_end  # this exits the global if statement
    else:  # find not exact C box (up to max 3 substitution allowed)
        for sub in range(1, int((len_c_box-1)/2 + 1)):  # iterate over 1 to 3 substitutions allowed
            c_motif = regex.findall("((A|G)UGAUGA){s<="+str(sub)+"}", first_20, overlapped=True)
            if len(c_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                c_motif = c_motif[0][0]  # if multiple C boxes found, keep the the C box closest to 5' end
                c_start = first_20.find(c_motif) + 1
                c_end = c_start + len(c_motif) - 1
                return c_motif, c_start, c_end  # this exits the global else statement
        # If no C box is found, return NNNNNNN and 0, 0 as C box sequence, start and end
        c_motif, c_start, c_end = 'NNNNNNN', 0, 0
        return c_motif, c_start, c_end


def find_d_box(seq):
    """ Find exact D box (CUGA), if not present, find D box with 1 or max 2
        substitutions. Return also the start and end position of that box as
        1-based values. If no D box is found, return a 'NNNN' empty D box and 0
        as start and end of D box."""
    first_20, last_20, length_seq = cut_sequence(seq)
    len_d_box = 4
    # First, find exact D box (CUGA) within 20 nt of the snoRNA 3' end
    if re.search('CUGA', last_20) is not None:  # find exact D box
        *_, last_possible_d = re.finditer('CUGA', last_20)
        d_motif = last_possible_d.group(0)  # if multiple exact D boxes found, keep the D box closest to 3' end
        d_start = (length_seq - 20) + last_possible_d.start() + 1
        d_end = (length_seq - 20) + last_possible_d.end()
        return d_motif, d_start, d_end
    else:  # find not exact D box (up to max 50% of substitution allowed (i.e. 2 nt))
        for sub in range(1, int(len_d_box/2 + 1)):  # iterate over 1 to 2 substitutions allowed
            d_motif = regex.findall("(CUGA){s<="+str(sub)+"}", last_20, overlapped=True)
            if len(d_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                d_motif = d_motif[-1]  # if multiple D boxes found, keep the the D box closest to 3' end
                d_start = (length_seq - 20) + last_20.rindex(d_motif) + 1
                d_end = d_start + len(d_motif) - 1
                return d_motif, d_start, d_end  # this exits the global else statement
        # If no D box is found, return NNNN and 0, 0 as D box sequence, start and end
        d_motif, d_start, d_end = 'NNNN', 0, 0
        return d_motif, d_start, d_end

def find_aca(line):
    """ Find the most downstream ACA motif in the last 10 nt of H/ACA box
        snoRNAs."""
    last_10 = line[-10:]
    length_seq = len(line)
    if re.search('ACA', last_10) is not None:  # find exact ACA box
        *_, last_possible_aca = re.finditer('ACA', last_10)
        aca_motif = last_possible_aca.group(0)  # if multiple exact ACA boxes found, keep the ACA box closest to 3' end
        aca_start = (length_seq - 10) + last_possible_aca.start() + 1  # 1-based position
        aca_end = (length_seq - 10) + last_possible_aca.end()  # 1-based
    else:  # if no ACA is found
        aca_motif, aca_start, aca_end = 'NNN', 0, 0

    return aca_motif, aca_start, aca_end



# Loading file and dfs
sno_fasta = snakemake.input.sno_fasta
tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
rna_central_df = pd.read_csv(snakemake.input.rna_central_df, sep='\t')
id_conversion_df = pd.read_csv(snakemake.input.id_conversion_df, sep='\t',
                        names=['rna_central_id', 'source', 'ensembl_transcript_id',
                                'taxid', 'biotype', 'ensembl_gene_id'])
mouse_gtf = pd.read_csv(snakemake.input.gtf, sep='\t', skiprows=5,
                        names=['chr', 'source', 'feature', 'start', 'end', 'dot',
                                'strand', 'dot2', 'attributes'])
mouse_gtf = mouse_gtf[mouse_gtf['feature'] == 'gene']

# Remove ".number" a the end of ensembl gene id
id_conversion_df[['ensembl_gene_id', 'suffix']] = id_conversion_df['ensembl_gene_id'].str.split('.', expand=True)

# Select only snoRNAs in the ensembl gtf
sno_only = mouse_gtf[mouse_gtf['attributes'].str.contains('gene_biotype "snoRNA"')]
sno_attributes = list(sno_only['attributes'])
sno_ids = [attribute.split(';') for attribute in sno_attributes]
sno_ids = [item for sublist in sno_ids for item in sublist]  # flatten list of list into one list
sno_ids = [attr.split(' "')[1].strip('"') for attr in sno_ids if 'gene_id' in attr]  # remove 'gene_id 'and '"'
sno_tpm_df = tpm_df[tpm_df['gene_id'].isin(sno_ids)]

# Dict of corresponding ensembl and rnacentral ids
id_mapping = dict(zip(id_conversion_df.ensembl_gene_id, id_conversion_df.rna_central_id))

# Find if sno_type is known for mouse snoRNAs in RNAcentral description
rna_central_df.loc[rna_central_df['description'].str.contains('C/D|SNORD|Snord|U3'), 'snoRNA_type'] = 'C/D'
rna_central_df.loc[rna_central_df['description'].str.contains('H/ACA|SNORA|Snora|ACA'), 'snoRNA_type'] = 'H/ACA'
sno_type_dict = dict(zip(rna_central_df.upi, rna_central_df.snoRNA_type))  # where upi is the column of RNA central id

# Create rna_central_id and snoRNA_type columns
sno_tpm_df['rna_central_id'] = sno_tpm_df['gene_id'].map(id_mapping)
sno_tpm_df['snoRNA_type'] = sno_tpm_df['rna_central_id'].map(sno_type_dict)

# For snoRNAs without snoRNA type from RNAcentral description, find if their gene_name contains this information directly
sno_tpm_df.loc[(sno_tpm_df['snoRNA_type'].isna()) & (sno_tpm_df['gene_name'].str.contains('Snord|SNORD|U8|U3')), 'snoRNA_type'] = 'C/D'
sno_tpm_df.loc[(sno_tpm_df['snoRNA_type'].isna()) & (sno_tpm_df['gene_name'].str.contains('Snora|SNORA')), 'snoRNA_type'] = 'H/ACA'

# Drop all snoRNAs where we couldn't find the snoRNA type (158 snoRNAs)
sno_tpm_df = sno_tpm_df.dropna(subset=['snoRNA_type'])

# For the snoRNAs without names telling if they are C/D or H/ACA, infer from their sequence if they have the canonical box motifs
'''
no_snoRNA_type_ids = list(sno_tpm_df[~sno_tpm_df['snoRNA_type'].isin(['C/D', 'H/ACA'])].gene_id)

# Create dictionary of ensembl_id: RNA sequence
seq_dict = {}
with open(sno_fasta, 'r') as f:
    temp_id = ''
    for line in f:
        if line.startswith('>'):
            id = line.strip('\n').strip('>')
            temp_id = id
        else:
            sno_sequence = line.strip('\n').replace('T', 'U')
            seq_dict[temp_id] = sno_sequence

# First, search for C and D boxes; if not present, search for ACA box
sno_type_search_dict = dict(zip(sno_tpm_df.gene_id, sno_tpm_df.snoRNA_type))
sno_type_search_dict = {k: sno_type_search_dict[k] for k in sno_type_search_dict if sno_type_search_dict[k] != 'NaN'}  # remove sno that have NaN as snoRNA type
other = 0
for sno_id in no_snoRNA_type_ids:
    snoRNA_seq = seq_dict[sno_id]
    c_motif, c_start, c_end = find_c_box(snoRNA_seq)
    d_motif, d_start, d_end = find_d_box(snoRNA_seq)
    if (c_motif == "NNNNNNN") | (d_motif == "NNNN"):  # if we don't find either a C or D box
        aca_motif, aca_start, aca_end = find_aca(snoRNA_seq)
        if aca_motif == "ACA":  # this is a ACA
            sno_type_search_dict[sno_id] = "H/ACA"
        else:
            other += 1
            print(f'{sno_id} is of unknown snoRNA type')
    else: # this is a C/D snoRNA
        sno_type_search_dict[sno_id] = "C/D"
print(f'{other} snoRNAs are of unknown snoRNA type. They are thereby excluded of downstream analyses.')

# Add the snoRNA type to these snoRNAs in sno_tpm_df and exclude the remaining snoRNAs that do not have a snoRNA type (only 25 snoRNAs for Mus musculus)
sno_tpm_df['snoRNA_type'] = sno_tpm_df['gene_id'].map(sno_type_search_dict)
sno_tpm_df = sno_tpm_df.dropna(subset=['snoRNA_type'])
len_sno_tpm_df = len(sno_tpm_df)
print(f'The snoRNA type (C/D or H/ACA) was found for {len_sno_tpm_df} snoRNAs. These snoRNAs are included in downstream analyses.')
'''
sno_tpm_df.to_csv(snakemake.output.snoRNA_type_df, sep='\t', index=False)
