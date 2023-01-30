#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

""" Determine the intersection between a bed file of all snoRNAs and a bedGraph
    of the phastCons conservation score for each nucleotide in the human genome.
    Then compute the average conservation score per snoRNA (over all the nt
    composing the snoRNA that have a conservation score associated)."""

col = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature',
        'dot2', 'gene_info']
sno_bed = pd.read_csv(snakemake.input.sno_bed, sep='\t', names=col)  # generated with gtf_to_bed

# Get bed of 100 nt upstream of snoRNA instead of snoRNA itself (account for strandness)
upstream = sno_bed.copy()
upstream['start_up'] = upstream['start'] - 1
upstream.loc[upstream['strand'] == '+', 'start'] = upstream['start_up'] - 100
upstream.loc[upstream['strand'] == '+', 'end'] = upstream['start_up']

upstream['end_up'] = upstream['end'] + 1
upstream.loc[upstream['strand'] == '-', 'end'] = upstream['end_up'] + 100
upstream.loc[upstream['strand'] == '-', 'start'] = upstream['end_up']
upstream = upstream[col]
upstream.to_csv('upstream_temp.bed', sep='\t', index=False, header=False)

def intersection(sno_bed, phastcons_bedgraph, output_path):
    """ Get the intersection between the snoRNA bed file and the sorted
        conservation bedgraph file."""
    a = BedTool(sno_bed)
    intersection = a.intersect(phastcons_bedgraph, wb=True, sorted=True).saveas(output_path)

def average_score(intersect_df_path, sno_bed_df, output_path):
    """Get the average conservation score per snoRNA."""
    intersect_df = pd.read_csv(intersect_df_path, sep='\t', names=['chr','start',
                                'end', 'gene_id', 'dot', 'strand', 'source', 'feature',
                                'dot2', 'gene_info', 'chr_bg', 'start_bg', 'end_bg', 'conservation_score'])

    # Groubpy every line (nucleotide) corresponding to a snoRNA and get the average conservation score across all these nt
    sno_nt = intersect_df.groupby(['gene_id'])['conservation_score'].mean()

    # Some snoRNAs are composed of nt that are not present in the conservation
    # bedgraph (only 2 snoRNAs don't have conservation scores at all; other snoRNAs might have a few nt missing a
    # conservation score, but these null values are not considered in the average
    # since they are not present in the intersect_df). For the 2 snoRNAs missing conservation score, we consider it as a score of 0
    sno_bed_df = sno_bed_df[['gene_id']]

    final_df = sno_bed_df.merge(sno_nt, how='left', left_on='gene_id', right_on='gene_id')
    final_df = final_df.fillna(0)
    final_df.to_csv(output_path, index=False, sep='\t')

def main(sno_bed_path, phastcons_bedgraph_path, output_path_intersection,
            sno_bed_df, output_path_final_df):
    intersect_df = intersection(sno_bed_path, phastcons_bedgraph_path, output_path_intersection)
    final_conservation_df = average_score(output_path_intersection, sno_bed_df, output_path_final_df)



# Get snoRNA average conservation across 100 vertebrates
main(snakemake.input.sno_bed, snakemake.input.phastcons_bg,
    snakemake.output.intersection_sno_conservation, sno_bed, snakemake.output.sno_conservation)

# Get average conservation of the 100 nt upstream of the snoRNAs (~promoter region)
main('upstream_temp.bed', snakemake.input.phastcons_bg,
    snakemake.output.intersection_upstream_sno_conservation, upstream, snakemake.output.upstream_sno_conservation)

sp.call('rm upstream_temp.bed', shell=True)
