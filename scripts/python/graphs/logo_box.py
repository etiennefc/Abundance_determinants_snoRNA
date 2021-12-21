#!/usr/bin/python3
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

""" Create logo of C and D boxes from fasta of either expressed or
    not expressed C/D box snoRNAs."""

fastas = snakemake.input.box_fastas
outputs = snakemake.output.box_logos

# Get all box sequences (not sno_id) in a list
for fasta in fastas:
    # Get the name of the box and if in expressed or not expressed snoRNAs to redirect figure to correct output
    ab_status_box = fasta.split('/')[-1].rstrip('.fa')
    output = [path for path in outputs if ab_status_box in path][0]
    with open(fasta) as f:
        raw_seqs = f.readlines()
    seqs = [seq.strip() for seq in raw_seqs if '>' not in seq]


    #Get a count and probability matrix to create the logo
    counts_matrix = logomaker.alignment_to_matrix(seqs)
    prob_matrix = logomaker.transform_matrix(counts_matrix, from_type='counts',
                                            to_type='probability')
    rc = {'ytick.labelsize': 32}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    logo = logomaker.Logo(prob_matrix, color_scheme='classic')
    logo.ax.set_ylabel("Frequency", fontsize=35)
    plt.savefig(output, bbox_inches='tight', dpi=600)
