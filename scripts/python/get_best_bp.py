#!/usr/bin/python3
import pandas as pd

""" Select only the predicted branch point nucleotide with the highest
    probability per HG intron and compute its distance to the snoRNA. The columns
    in the output df are seqnames (chr), start (start position of the window,
    i.e always 44 nt upstream of the 3'exon), end (end position of the window,
    i.e always 18 nt upstream of the 3'exon), strand, gene_id, transcript_id,
    exon_3prime (exon_id of 3' exon), exon_5prime (exon_id of 5' exon),
    exon_number (exon number of exon downstream of the branch point),
    intron_number (intron in which the branch point is located i.e the number of
    the upstream exon), test_site (position of the branch_point), seq_pos0 (nt
    at the branch point), branchpoint_prob (branchpoint probability computed by
    branchpointer), U2_binding_energy (U2 binding energy to the branch point
    computed by branchpointer) and bp_to_3prime (distance of branch point to
    3' exon)."""

bp_total_df = pd.read_csv(snakemake.input.bp_distance_total_df)
sno_location_df = pd.read_csv(snakemake.input.sno_location_df, sep='\t')
sno_overlap_df = pd.read_csv(snakemake.params.sno_overlap_df, sep='\t')

#Split bp_total_df into a SNHG14 df and a 'all other host genes' (HG) df
snhg14_bp = bp_total_df[bp_total_df['transcript_id'] == "NR_146177.1"]
all_other_hg_bp = bp_total_df[bp_total_df['transcript_id'] != "NR_146177.1"]

## Grouby df by chr, exon3prime and exon5prime of window and take only the line with the highest probability of branch point per group
#For all HG except SNHG14
all_other_hg_bp['max'] = all_other_hg_bp.groupby(['seqnames', 'exon_3prime', 'exon_5prime'])['branchpoint_prob'].transform('max')
all_other_hg_bp = all_other_hg_bp[all_other_hg_bp['max'] == all_other_hg_bp['branchpoint_prob']]  # Select only the max probability branch point

#For SNHG14, groupby exon number only (since no exon ids are available in the refseq gtf)
snhg14_bp['max'] = snhg14_bp.groupby('exon_number')['branchpoint_prob'].transform('max')
snhg14_bp = snhg14_bp[snhg14_bp['max'] == snhg14_bp['branchpoint_prob']]  # Select only the max probability branch point


# Calculate bp_to_3prime distance and intron number for both dfs and concat the resulting dfs into one df
all_other_hg_bp['bp_to_3prime'] = all_other_hg_bp['end'] + 18 - all_other_hg_bp['test_site']
all_other_hg_bp['intron_number'] = all_other_hg_bp['exon_number'] - 1

all_other_hg_bp = all_other_hg_bp[['seqnames', 'start', 'end', 'strand', 'gene_id',
                'transcript_id', 'exon_3prime', 'exon_5prime',
                'exon_number', 'intron_number', 'test_site', 'seq_pos0',
                'branchpoint_prob', 'U2_binding_energy', 'bp_to_3prime']]

snhg14_bp['bp_to_3prime'] = snhg14_bp['end'] + 18 - snhg14_bp['test_site']
snhg14_bp['intron_number'] = snhg14_bp['exon_number'] - 1

snhg14_bp = snhg14_bp[['seqnames', 'start', 'end', 'strand', 'gene_id',
                'transcript_id', 'exon_3prime', 'exon_5prime',
                'exon_number', 'intron_number', 'test_site', 'seq_pos0',
                'branchpoint_prob', 'U2_binding_energy', 'bp_to_3prime']]

bp_distance_simple = pd.concat([all_other_hg_bp, snhg14_bp], ignore_index=True)
bp_distance_simple.to_csv(snakemake.output.bp_distance_simple, index=False, sep='\t')

# Merge bp_distance_simple df to sno_location_df and calculate the distance between intronic snoRNAs and the predicted branch_point in their intron
sno_location_df['intron_number'] = sno_location_df['intron_number'].astype(int)
sno_bp_df = sno_location_df.merge(bp_distance_simple[['transcript_id', 'intron_number', 'bp_to_3prime']],
                how='left', left_on=['transcript_id_host', 'intron_number'], right_on=['transcript_id', 'intron_number'])

sno_bp_df['dist_to_bp'] = sno_bp_df['distance_downstream_exon'] - sno_bp_df['bp_to_3prime']

# Replace bp_to_3prime and dist_to_bp by 0 for snoRNAs that overlap a HG exon
overlapping_sno = list(sno_overlap_df[~sno_overlap_df['hg_overlap'].isna()].sno)
overlapping_sno.append('NR_000026') #this snoRNA is at 2 nt from its downstream exon, so the dist_to_bp would be negative; so we consider it as overlaping also
for i, sno_id in enumerate(overlapping_sno):
    sno_bp_df.loc[sno_bp_df.gene_id_sno == sno_id, 'bp_to_3prime'] = 0
    sno_bp_df.loc[sno_bp_df.gene_id_sno == sno_id, 'dist_to_bp'] = 0

sno_bp_df.drop(['transcript_id_x', 'transcript_id_y'], axis=1, inplace=True)

sno_bp_df.to_csv(snakemake.output.sno_distance_bp, index=False, sep='\t')
