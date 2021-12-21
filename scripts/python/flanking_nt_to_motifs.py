#!/usr/bin/python3
import pandas as pd

expressed_cd_fa = snakemake.input.expressed_cd_fa
not_expressed_cd_fa = snakemake.input.not_expressed_cd_fa
expressed_haca_fa = snakemake.input.expressed_haca_fa
not_expressed_haca_fa = snakemake.input.not_expressed_haca_fa

c_d_box_expressed = pd.read_csv(snakemake.input.c_d_box_location_expressed, sep='\t')
c_d_box_not_expressed = pd.read_csv(snakemake.input.c_d_box_location_not_expressed, sep='\t')
h_aca_box_expressed = pd.read_csv(snakemake.input.h_aca_box_location_expressed, sep='\t')
h_aca_box_not_expressed = pd.read_csv(snakemake.input.h_aca_box_location_not_expressed, sep='\t')

def dict_to_fasta(dictio, output_path):
    """ Use dict to create fasta (keys being used for ID lines and
        values being sequence lines)."""
    with open(output_path, 'w') as f:
        for sno_id, sequence in dictio.items():
            f.write(f'>{sno_id}\n')
            f.write(f'{sequence}\n')



def flanking_nt_to_c_box(position_df, input_fasta, output_fasta, box_type):
    """ Find start and end of C box in position_df and extract 3nt flanking up- and downstream
        of C box. Return flanking nt + C box sequence as fasta. If C box is directly at the 5'
        start of the snoRNA, return 'nnn' as the flanking nt. Also, return flanking nt in
        lowercase and motif as uppercase. This function can also be applied to C' and H boxes
        using box_type."""
    len_motif_dict = {"C": "nnnNNNNNNNnnn", "C_prime": "nnnNNNNNNNnnn", "H": "nnnNNNNNNnnn"}
    flanking_motif_dict = {}
    position_df = position_df.set_index('gene_id')
    start_end_dict = position_df[[f'{box_type}_start', f'{box_type}_end']].to_dict('index')  # these positions are 1-based
    with open(input_fasta, 'r') as f:
        sno_id = ''
        for line in f:
            if line.startswith('>'):
                id = line.lstrip('>').rstrip('\n')
                sno_id = id
            else:
                seq = line.rstrip('\n')
                start, end = start_end_dict[sno_id][f'{box_type}_start'], start_end_dict[sno_id][f'{box_type}_end']
                if (start == 0) & (end == 0):  # i.e. no C box was found
                    flanking_motif_dict[sno_id] = len_motif_dict[box_type]
                elif start >= 4:  # i.e. a C box was found after the first 3 nt
                    flanking_motif = seq[start-4:end+3]  # -4 and +3 to convert to 0-based indexing
                    left, right = flanking_motif[0:3].lower(), flanking_motif[-3:].lower()
                    flanking_motif_dict[sno_id] = left + flanking_motif[3:-3] + right
                else:  # i.e. a C box was found but included within the first 3 nt
                    n = {0: 'nnn', 1: 'nn', 2:'n'}
                    flanking_motif = seq[0:end+3]
                    flanking_motif = n[start-1] + flanking_motif
                    left, right = flanking_motif[0:3].lower(), flanking_motif[-3:].lower()
                    flanking_motif_dict[sno_id] = left + flanking_motif[3:-3] + right

    dict_to_fasta(flanking_motif_dict, output_fasta)


def flanking_nt_to_d_box(position_df, input_fasta, output_fasta, box_type):
    """ Find start and end of D box in position_df and extract 3nt flanking up- and downstream
        of D box. Return flanking nt + D box sequence as fasta. If D box is directly at the 3'
        start of the snoRNA, return 'nnn' as the flanking nt. Also, return flanking nt in
        lowercase and motif as uppercase. This function can also be applied to D' and ACA boxes
        using box_type."""
    len_motif_dict = {"D": "nnnNNNNnnn", "D_prime": "nnnNNNNnnn", "ACA": "nnnNNNnnn"}
    flanking_motif_dict = {}
    position_df = position_df.set_index('gene_id')
    start_end_dict = position_df[[f'{box_type}_start', f'{box_type}_end']].to_dict('index')  # these positions are 1-based
    with open(input_fasta, 'r') as f:
        sno_id = ''
        for line in f:
            if line.startswith('>'):
                id = line.lstrip('>').rstrip('\n')
                sno_id = id
            else:
                seq = line.rstrip('\n')
                start, end = start_end_dict[sno_id][f'{box_type}_start'], start_end_dict[sno_id][f'{box_type}_end']
                if (start == 0) & (end == 0):  # i.e. no D box was found
                    flanking_motif_dict[sno_id] = len_motif_dict[box_type]
                elif end <= len(seq) - 3:  # i.e. a D box was found before the last 3 nt
                    flanking_motif = seq[start-4:end+3]  # -4 and +3 to convert to 0-based indexing
                    left, right = flanking_motif[0:3].lower(), flanking_motif[-3:].lower()
                    flanking_motif_dict[sno_id] = left + flanking_motif[3:-3] + right
                else:  # i.e. a D box was found but included within the last 3 nt
                    n = {0: 'nnn', 1: 'nn', 2:'n'}
                    flanking_motif = seq[start-4:]
                    flanking_motif = flanking_motif + n[len(seq) - end]
                    left, right = flanking_motif[0:3].lower(), flanking_motif[-3:].lower()
                    flanking_motif_dict[sno_id] = left + flanking_motif[3:-3] + right

    dict_to_fasta(flanking_motif_dict, output_fasta)


# Create fasta of motif and flanking nt for expressed and not expressed C/D box snoRNAs
flanking_nt_to_c_box(c_d_box_expressed, expressed_cd_fa, snakemake.output.c_expressed, 'C')
flanking_nt_to_c_box(c_d_box_expressed, expressed_cd_fa, snakemake.output.c_prime_expressed, 'C_prime')
flanking_nt_to_d_box(c_d_box_expressed, expressed_cd_fa, snakemake.output.d_expressed, 'D')
flanking_nt_to_d_box(c_d_box_expressed, expressed_cd_fa, snakemake.output.d_prime_expressed, 'D_prime')

flanking_nt_to_c_box(c_d_box_not_expressed, not_expressed_cd_fa, snakemake.output.c_not_expressed, 'C')
flanking_nt_to_c_box(c_d_box_not_expressed, not_expressed_cd_fa, snakemake.output.c_prime_not_expressed, 'C_prime')
flanking_nt_to_d_box(c_d_box_not_expressed, not_expressed_cd_fa, snakemake.output.d_not_expressed, 'D')
flanking_nt_to_d_box(c_d_box_not_expressed, not_expressed_cd_fa, snakemake.output.d_prime_not_expressed, 'D_prime')

# Create fasta of motif and flanking nt for expressed and not expressed H/ACA box snoRNAs
flanking_nt_to_c_box(h_aca_box_expressed, expressed_haca_fa, snakemake.output.h_expressed, 'H')
flanking_nt_to_d_box(h_aca_box_expressed, expressed_haca_fa, snakemake.output.aca_expressed, 'ACA')

flanking_nt_to_c_box(h_aca_box_not_expressed, not_expressed_haca_fa, snakemake.output.h_not_expressed, 'H')
flanking_nt_to_d_box(h_aca_box_not_expressed, not_expressed_haca_fa, snakemake.output.aca_not_expressed, 'ACA')
