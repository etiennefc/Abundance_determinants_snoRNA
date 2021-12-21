#!/usr/bin/python3
import pandas as pd
import re
import regex
from functools import reduce

""" Find H and ACA boxes of each snoRNA (if they exist) and their position. The
    found boxes need to be exact (no substitution allowed) given the already not
    so precise H and ACA motifs."""
dot_bracket = snakemake.input.dot_bracket
expressed_haca = snakemake.input.expressed_haca
not_expressed_haca = snakemake.input.not_expressed_haca
output_expressed_haca = snakemake.output.h_aca_box_location_expressed
output_not_expressed_haca = snakemake.output.h_aca_box_location_not_expressed

# Get dot bracket for all snoRNAs inside dict
structure = {}
with open(dot_bracket, 'r') as f:
    sno_id = ''
    for line in f:
        line = line.strip('\n')
        if line.startswith('>'):
            id = line.strip('>')
            sno_id = id
        elif line.startswith(('(', ')', '.')):
            dot_bracket_structure = line.split(' ')[0]
            structure[sno_id] = dot_bracket_structure


def generate_df(dictio, motif_name):
    """ From a dictionary containing the motif and start/end per H/ACA id"""
    # Create dataframe from box_dict
    box = pd.DataFrame.from_dict(dictio, orient='index')
    box.columns = [f'{motif_name}_sequence', f'{motif_name}_start', f'{motif_name}_end']
    box = box.reset_index()
    box = box.rename(columns={"index": "gene_id"})
    return box


def find_h_box(fasta, dot_bracket_dict):
    """ Find potential Hinge region(s) (if it exists) in dot bracket of H/ACA
        snoRNAs. This region in the dot bracket is represented by ')....('. If
        an ANANNA is present, return that box; otherwise return NNNNNN (no
        mismatch allowed since H box is already not so well defined)."""
    h_box_dict = {}
    with open(fasta, 'r') as f:
        haca_id = ''
        for line in f:
            line = line.strip('\n')
            if line.startswith('>'):
                id = line.strip('>')
                haca_id = id
            else:
                dot_bracket_temp = dot_bracket_dict[haca_id]
                if re.search('\)\.{6,}\(', dot_bracket_temp) is not None:  # find possible hinge regions (at least 6 unpaired nucleotides in the middle of the snoRNA)
                    hinges = re.finditer('\)\.{6,}\(', dot_bracket_temp)
                    h_dot_bracket = re.findall('\)\.{6,}\(', dot_bracket_temp)
                    start = [h.start(0) for h in hinges]
                    end = [start[i] + len(h) for i, h in enumerate(h_dot_bracket)]
                    for i, hinge in enumerate(h_dot_bracket):
                        seq = line[start[i]+1:end[i]-1]  # get hinge region unpaired nucleotides (only the '.', not the surrounding ')' or '(')
                        if re.search('A.A..A', seq) is not None:  # find the first (closest to 5') exact H box within possible hinge regions
                            h_motifs = re.findall('A.A..A', seq)
                            substart = seq.index(h_motifs[0])  # start of H box within the extracted hinge region
                            h_start = start[i] + 2 + substart  # +2 because +1 to be the first unpaired nt in the hinge region and +1 to be 1-based
                            h_box_dict[haca_id] = {'h_motif': h_motifs[0], 'h_start': h_start, 'h_end': h_start + 5}  # +5 nt after the first a in H box
                            break
                    if haca_id not in h_box_dict.keys():  # if no H box was found within the found hinge region(s)
                        h_motif, h_start, h_end = 'NNNNNN', 0, 0
                        h_box_dict[haca_id] = {'h_motif': h_motif, 'h_start': h_start, 'h_end': h_end}

                else:  # if no unpaired hinge region was found, return NNNNNN H box and 0 as start and end
                    h_motif, h_start, h_end = 'NNNNNN', 0, 0
                    h_box_dict[haca_id] = {'h_motif': h_motif, 'h_start': h_start, 'h_end': h_end}

    h_box_df = generate_df(h_box_dict, 'H')
    return h_box_df




def find_aca(fasta):
    """ Find the most downstream ACA motif in the last 10 nt of H/ACA box
        snoRNAs."""
    aca_box_dict = {}
    with open(fasta, 'r') as f:
        haca_id = ''
        for line in f:
            line = line.strip('\n')
            if line.startswith('>'):
                id = line.strip('>')
                haca_id = id
            else:
                last_10 = line[-10:]
                length_seq = len(line)
                if re.search('ACA', last_10) is not None:  # find exact ACA box
                    *_, last_possible_aca = re.finditer('ACA', last_10)
                    aca_motif = last_possible_aca.group(0)  # if multiple exact ACA boxes found, keep the ACA box closest to 3' end
                    aca_start = (length_seq - 10) + last_possible_aca.start() + 1  # 1-based position
                    aca_end = (length_seq - 10) + last_possible_aca.end()  # 1-based
                    aca_box_dict[haca_id] = {'aca_motif': aca_motif, 'aca_start': aca_start, 'aca_end': aca_end}
                else:  # if no ACA is found
                    aca_box_dict[haca_id] = {'aca_motif': 'NNN', 'aca_start': 0, 'aca_end': 0}

    aca_box_df = generate_df(aca_box_dict, 'ACA')
    return aca_box_df


def find_all_boxes(fasta, dot_bracket_dict, path):
    """ Find H and ACA boxes in given fasta using find_h_box and find_aca and
        concat resulting dfs horizontally."""
    df_h = find_h_box(fasta, dot_bracket_dict)
    df_aca = find_aca(fasta)

    df_final = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'],
                                            how='outer'),
                                            [df_h, df_aca])
    df_final.to_csv(path, index=False, sep='\t')



find_all_boxes(expressed_haca, structure, output_expressed_haca)
find_all_boxes(not_expressed_haca, structure, output_not_expressed_haca)
