#!/usr/bin/python3
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from math import log2
from scipy.stats import entropy, kstest

""" Create logo of C and D boxes from fasta of either expressed or
    not expressed C/D box snoRNAs."""

fastas = snakemake.input.box_fastas
logo_outputs = snakemake.output.box_logos
pie_outputs = snakemake.output.pie_logos
color_dict = snakemake.params.color_dict

def get_entropy(proba_matrix):
    """ Compute the entropy of a given logo from a proba_matrix (each column is
        a nucleotide (A, U, C or G), each line is the position in the logo)."""
    cumulative_entropy = []
    for i, row in proba_matrix.iterrows():
        vals = list(row)
        entropy_per_position = entropy(vals, base=2)
        cumulative_entropy.append(entropy_per_position)
    print(sum(cumulative_entropy))
    return sum(cumulative_entropy)

def make_autopct(values):
    """ Create function to return % in pie chart"""
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d} \n ({p:.1f}%)'.format(p=pct,v=val)
    return my_autopct

def pie_simple_annot(count_list, colors, annotation, path, **kwargs):
    """
    Creates a pie chart from a simple list of values and add an annotation
    (ex: an entropy value) next to the graph.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))

    ax.pie(count_list, colors=list(colors.values()), pctdistance=0.5,
           textprops={'fontsize': 35}, autopct=make_autopct(count_list), **kwargs)
    fig.suptitle(annotation, y=0.9, fontsize=18)
    plt.legend(labels=colors.keys(), loc='upper right',
                bbox_to_anchor=(1, 1.18), prop={'size': 35})
    plt.savefig(path, dpi=600)


# Get all box sequences (not sno_id) in a list
f_obs_dict = {}
for fasta in fastas:
    # Get the name of the box and if in expressed or not expressed snoRNAs to redirect figure to correct output
    ab_status_box = fasta.split('/')[-1].rstrip('.fa')
    logo_output = [path for path in logo_outputs if ab_status_box in path][0]
    pie_output = [path for path in pie_outputs if ab_status_box in path][0]
    with open(fasta) as f:
        raw_seqs = f.readlines()
    seqs = [seq.strip() for seq in raw_seqs if '>' not in seq]
    seqs_wo_blank = [seq for seq in seqs if 'NNN' not in seq]  # remove blank (NNNN) sequences
    len_seqs, len_seqs_wo_blank = len(seqs), len(seqs_wo_blank)

    #Get a count and probability matrix to create the logo
    counts_matrix = logomaker.alignment_to_matrix(seqs_wo_blank)
    prob_matrix = logomaker.transform_matrix(counts_matrix, from_type='counts',
                                            to_type='probability')

   # Create logo wo blanks
    rc = {'ytick.labelsize': 32}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    logo = logomaker.Logo(prob_matrix, color_scheme='classic')
    logo.ax.set_ylabel("Frequency", fontsize=35)
    plt.savefig(logo_output, bbox_inches='tight', dpi=600)

    # Compute entropy of each logo (blanks removed)
    entropy_ = str(get_entropy(prob_matrix))

    print(ab_status_box, prob_matrix)

    # Get the observed frequency into a flattened list
    l = prob_matrix.values.tolist()
    f_obs = [j for sublist in l for j in sublist]
    f_obs_dict[ab_status_box] = f_obs
    
    # Create pie chart of found vs not found motif per box
    percent = [(len_seqs_wo_blank/len_seqs)*100, ((len_seqs - len_seqs_wo_blank)/len_seqs)*100]
    pie_simple_annot(percent, color_dict, f'{entropy_} bits', pie_output)

# Compute Kolmogorov-Smirnov test to see if the statistical significance of box degeneration between expressed/not expressed snoRNAs
res = kstest(f_obs_dict['not_expressed_c_box'], f_obs_dict['expressed_c_box'])
print(res)
res = kstest(f_obs_dict['not_expressed_d_box'], f_obs_dict['expressed_d_box'])
print(res)
res = kstest(f_obs_dict['not_expressed_c_prime_box'], f_obs_dict['expressed_c_prime_box'])
print(res)
res = kstest(f_obs_dict['not_expressed_d_prime_box'], f_obs_dict['expressed_d_prime_box'])
print(res)
res = kstest(f_obs_dict['not_expressed_aca_box'], f_obs_dict['expressed_aca_box'])
print(res)
res = kstest(f_obs_dict['not_expressed_h_box'], f_obs_dict['expressed_h_box'])
print(res)

