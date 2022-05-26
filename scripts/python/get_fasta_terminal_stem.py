#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool

""" Get the sequences composing the potential terminal stems of snoRNAs into a
    combined fasta file. Each entry correspond to the left flanking region in
    reverse order followed by a '&' and the right flanking region in reverse
    order. This orders ensures that RNAcofold respects the 5'->3' opposite base
    pairing composing the potential terminal stems and create representative
    structure graphs. On the resulting graphs, the green sequence corresponds to
    the left flanking region and the red sequence corresponds to the right
    flanking region (the 5' of these sequences being at the extremity with no
    full dot inside the nucleotide on the graph; the 3' of these sequences being
    at the extermitiy where there is a full dot inside the nucleotide). """
col_names = ['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'source', 'feature', 'score2', 'characteristics']
cd_left = BedTool(snakemake.input.flanking_cd_left)
cd_right = BedTool(snakemake.input.flanking_cd_right)
haca_left = BedTool(snakemake.input.flanking_haca_left)
haca_right = BedTool(snakemake.input.flanking_haca_right)

# Verify if there is a left and right sequence for all snoRNAs (sometimes the
# snoRNA is at the end of a chr, which makes it impossible to find a right sequence)
cd_left_df = pd.read_table(cd_left.fn, names=col_names)
cd_right_df = pd.read_table(cd_right.fn, names=col_names)
if list(cd_left_df.gene_id) != list(cd_right_df.gene_id):
    diff_id = list(set(list(cd_left_df.gene_id)) - set(list(cd_right_df.gene_id)))
    diff_id2 = list(set(list(cd_right_df.gene_id)) - set(list(cd_left_df.gene_id)))
    diff_ids = diff_id + diff_id2
    print(diff_ids)
    cd_left_df = cd_left_df[~cd_left_df['gene_id'].isin(diff_ids)]
    cd_right_df = cd_right_df[~cd_right_df['gene_id'].isin(diff_ids)]
    cd_left = BedTool.from_dataframe(cd_left_df)
    cd_right = BedTool.from_dataframe(cd_right_df)

# Get the sequences of the extended flanking regions of C/D snoRNAs
fasta_cd_left = cd_left.sequence(fi=snakemake.input.genome_fasta, s=True)
fasta_cd_right = cd_right.sequence(fi=snakemake.input.genome_fasta, s=True)

# For C/D snoRNAs, open both left and right flanking regions fasta file and append the right region
# to the right of the left region as reverse order strings separated by a '&'; this combined sequence is used by RNAcofold.
# The attribute "seqfn" points to the new fasta file created by sequence().
seqs_cd = []
with open(fasta_cd_left.seqfn, 'r') as file_left, open(fasta_cd_right.seqfn, 'r') as file_right:
    for line_left, line_right in zip(file_left, file_right):
        line_left = str(line_left)
        line_right = str(line_right)
        if (not line_left.startswith('>')) & (not line_right.startswith('>')):
            # Reverse the order of both flanking sequences with [::-1]
            co_seq = line_left[::-1] + "&" + line_right[::-1]  # this order is how RNAcofold will accurately try to see the
                                                # best base pairing between the left and right extended flanking regions of snoRNAs
            co_seq = co_seq.replace('T', 'U')  # convert DNA to RNA
            co_seq = co_seq.replace('\n', '')  # remove new lines from string
            seqs_cd.append(co_seq)
print(len(seqs_cd))
# Create a dictionary of C/D snoRNAs id as keys and co_seq as values
cd_dictio = {}
cd_ids = pd.read_table(cd_left.fn, names=col_names)
for i, gene_id in enumerate(cd_ids['gene_id']):
    cd_dictio[gene_id] = seqs_cd[i]
cd_dictio = {'>'+ k: v for k, v in cd_dictio.items()}  # Add '>' in front of all sno id

# Append these C/D snoRNAs ids and co_seq into the output file
with open(snakemake.output.sequences, "a+") as file:  # a+ for append in new file
    for k, v in cd_dictio.items():
        file.write(k+'\n'+v+'\n')



# Get the sequences of the extended flanking regions of H/ACA snoRNAs
fasta_haca_left = haca_left.sequence(fi=snakemake.input.genome_fasta, s=True)
fasta_haca_right = haca_right.sequence(fi=snakemake.input.genome_fasta, s=True)

# For H/ACA snoRNAs, open both left and right flanking regions fasta file and append the right region
# to the right of the left region as reverse order strings separated by a '&'; this combined sequence is used by RNAcofold.
# The attribute "seqfn" points to the new fasta file created by sequence().
seqs_haca = []
with open(fasta_haca_left.seqfn, 'r') as file_left, open(fasta_haca_right.seqfn, 'r') as file_right:
    for line_left, line_right in zip(file_left, file_right):
        line_left = str(line_left)
        line_right = str(line_right)
        if (not line_left.startswith('>')) & (not line_right.startswith('>')):
            # Reverse the order of both flanking sequences with [::-1]
            co_seq = line_left[::-1] + "&" + line_right[::-1]  # this order is how RNAcofold will accurately try to see the
                                                # best base pairing between the left and right extended flanking regions of snoRNAs
            co_seq = co_seq.replace('T', 'U')  # convert DNA to RNA
            co_seq = co_seq.replace('\n', '')  # remove new lines from string
            seqs_haca.append(co_seq)

# Create a dictionary of H/ACA snoRNAs id as keys and co_seq as values
haca_dictio = {}
haca_ids = pd.read_table(haca_left.fn, names=col_names)
for i, gene_id in enumerate(haca_ids['gene_id']):
    haca_dictio[gene_id] = seqs_haca[i]
haca_dictio = {'>'+ k: v for k, v in haca_dictio.items()}  # Add '>' in front of all sno id

# Append these H/ACA snoRNAs ids and co_seq into the output file
with open(snakemake.output.sequences, "a+") as file:  # a+ for append in output file already containing C/D snoRNAs co_seq
    for k, v in haca_dictio.items():
        file.write(k+'\n'+v+'\n')
