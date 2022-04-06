#!/usr/bin/python3
import pandas as pd

sno_fasta = snakemake.input.sno_fasta
snodb = pd.read_csv(snakemake.input.snodb, sep='\t')
cd_output, haca_output = snakemake.output.cd_fasta, snakemake.output.haca_fasta

# Select either C/D or H/ACA box snoRNAs in fasta of all snoRNA sequences
cd_ids = snodb[snodb['sno_type'] == 'C/D']['gene_id_sno'].to_list()
haca_ids = snodb[snodb['sno_type'] == 'H/ACA']['gene_id_sno'].to_list()
cd_dict, haca_dict = {}, {}
with open(sno_fasta, 'r') as f:
    sno_id = ''
    for line in f:
        if line.startswith('>'):
                id = line.lstrip('>').rstrip('\n')
                sno_id = id
        else:
            seq = line.rstrip('\n')
            if sno_id in cd_ids:
                cd_dict[sno_id] = seq
            elif sno_id in haca_ids:
                haca_dict[sno_id] = seq

# Create fasta of C/D snoRNAs
with open(cd_output, 'w') as f:
    for sno_id, sequence in cd_dict.items():
        f.write(f'>{sno_id}\n')
        f.write(f'{sequence}\n')

# Create fasta of H/ACA snoRNAs
with open(haca_output, 'w') as f:
    for sno_id, sequence in haca_dict.items():
        f.write(f'>{sno_id}\n')
        f.write(f'{sequence}\n')



def generate_df_prime(fasta):
    """ From a fasta of snoRNA sequences, find a given motif (C' or D') using
        predefined function find_c_prime_d_prime_hamming and output the motif sequence,
        start and end as a df."""
    # Get motif, start and end position inside dict
    box_dict = {}
    with open(fasta, 'r') as f:
        sno_id = ''
        for line in f:
            if line.startswith('>'):
                id = line.lstrip('>').rstrip('\n')
                sno_id = id
            else:
                seq = line.rstrip('\n')
                c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end = find_c_prime_d_prime_hamming(seq)
                box_dict[sno_id] = [c_prime_motif, c_prime_start,
                                    c_prime_end, d_prime_motif, d_prime_start,
                                    d_prime_end]

    # Create dataframe from box_dict
    box = pd.DataFrame.from_dict(box_dict, orient='index',
                                columns=['C_prime_sequence', 'C_prime_start',
                                        'C_prime_end', 'D_prime_sequence',
                                        'D_prime_start', 'D_prime_end'])
    box = box.reset_index()
    box = box.rename(columns={"index": "gene_id"})
    return box
