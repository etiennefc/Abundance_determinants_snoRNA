#!/usr/bin/python3
import pandas as pd

c_d_box = pd.read_csv(snakemake.input.c_d_box_location, sep='\t')
h_aca_box = pd.read_csv(snakemake.input.h_aca_box_location, sep='\t')


def cols_to_dict(df, col1, col2):
    """ Convert two columns of df into dictionary (col1 as keys and col2 as values). """
    df = df[[col1, col2]]
    df = df.set_index(col1)
    dictio = df.to_dict('index')
    return dictio



def hamming(found_motif, consensus_motif):
    """ Find the hamming distance of a found motif compared to the consensus
        motif. """
    hamming = 0
    for i, char in enumerate(found_motif):
        if char != consensus_motif[i]:
            hamming += 1
    return hamming


def convert_h_box(h_motif):
    """ Convert H motif for hamming distance purposes. Convert NNNNNN into ZZZZZZ 
        (so that it is counted as totally different from ANANNA), and convert the 
        2nd, 4th and 5th nucleotide into a N if a H motif was found (so the 3 N 
        in the found motif are counted as matching with the consensus motif)"""
    if h_motif == 'NNNNNN':
        h_motif = 'ZZZZZZ'
    else:
        h_motif = h_motif[0] + 'N' + h_motif[2] + 'NN' + h_motif[5]
    return h_motif


# Create one dictionary per motif (where key/val are sno_id/motif_sequence)
motifs = ['D_sequence', 'C_sequence', 'D_prime_sequence', 'C_prime_sequence', 'H_sequence', 'ACA_sequence']
motif_dicts = []
for motif in motifs:
    if motif.startswith(('C','D')):
        d = cols_to_dict(c_d_box, 'gene_id', motif)
    else:
        d = cols_to_dict(h_aca_box, 'gene_id', motif)
    motif_dicts.append(d)


# Compute hamming distance for each motif
seq = {'D_sequence': 'CUGA', 'C_sequence': 'RUGAUGA', 'D_prime_sequence': 'CUGA', 
    'C_prime_sequence': 'RUGAUGA', 'H_sequence': 'ANANNA', 'ACA_sequence': 'ACA'}
hamming_dict = {}
for dictio in motif_dicts:  # iterate through every motif dictionary (one for C, one for D, one for C_prime, etc.)
    for sno_id, motif_dict in dictio.items():
        if sno_id not in hamming_dict.keys():  # integrates for the first time in hamming_dict a sno_id with a given motif hamming distance
            for k, motif in motif_dict.items():
                if k.startswith('C'):  # for C and C_prime motifs
                    if motif.startswith(('A', 'G')):  # convert first nt of C or C prime motif to R if it's A|G
                        motif = 'R' + motif[1:]
                        hamming_val = hamming(motif, seq[k])
                        hamming_dict[sno_id] = {k: hamming_val}
                    else:  # do not convert to R the first nt because it is not a A|G
                        hamming_val = hamming(motif, seq[k])
                        hamming_dict[sno_id] = {k: hamming_val}
                elif k.startswith('H'):  # convert H motif according to convert_h_box
                    motif = convert_h_box(motif)
                    hamming_val = hamming(motif, seq[k])
                    hamming_dict[sno_id] = {k: hamming_val}
                else:  # for D, D prime and ACA motifs
                    hamming_val = hamming(motif, seq[k])
                    hamming_dict[sno_id] = {k: hamming_val}
        else:  # integrates in hamming_dict the hamming distance of new motifs for sno_ids that are already present in the hamming_dict because of the precedent global if statement
            for k, motif in motif_dict.items():
                if k.startswith('C'):
                    if motif.startswith(('A', 'G')):
                        motif = 'R' + motif[1:]
                        hamming_val = hamming(motif, seq[k])
                        hamming_dict[sno_id][k] = hamming_val
                    else:
                        hamming_val = hamming(motif, seq[k])
                        hamming_dict[sno_id][k] = hamming_val
                elif k.startswith('H'):
                    motif = convert_h_box(motif)
                    hamming_val = hamming(motif, seq[k])
                    hamming_dict[sno_id][k] = hamming_val
                else:
                    hamming_val = hamming(motif, seq[k])
                    hamming_dict[sno_id][k] = hamming_val

# Create df, reorder cols, change col names and compute combined hamming distance (sum of hamming distance for all boxes of a given snoRNA)
df = pd.DataFrame.from_dict(hamming_dict, orient='index')
df = df.reset_index()
df = df.rename(columns={'index': 'gene_id'})
df = df[['gene_id', 'C_sequence', 'D_sequence', 'C_prime_sequence', 'D_prime_sequence', 'H_sequence', 'ACA_sequence']]
df.columns = ['gene_id', 'C_hamming', 'D_hamming', 'C_prime_hamming', 'D_prime_hamming', 'H_hamming', 'ACA_hamming']
df['C_D_and_prime_hamming'] = df['C_hamming'] + df['D_hamming'] + df['C_prime_hamming'] + df['D_prime_hamming'] 
df['H_ACA_hamming'] = df['H_hamming'] + df['ACA_hamming']
df['combined_box_hamming'] = df['C_D_and_prime_hamming'].fillna(df['H_ACA_hamming'])
df = df.drop(columns=['C_D_and_prime_hamming', 'H_ACA_hamming'])

df.to_csv(snakemake.output.hamming_distance_box_df, index=False, sep='\t')


