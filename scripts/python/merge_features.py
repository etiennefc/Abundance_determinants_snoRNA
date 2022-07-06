#!/usr/bin/python3
import pandas as pd
import numpy as np

""" Generate a merged dataframe of all snoRNA features and labels that will be
    used by the predictor. The labels are 'expressed' or 'not_expressed' using
    the column 'abundance_cutoff_2'. The features are: snoRNA length, type,
    target, host gene (HG) biotype and function, HG NMD susceptibility, HG
    promoter type, intron rank (from 5', 3' and relative to the total number of
    intron) and length in which the snoRNA is encoded, total number of intron of
    the snoRNA's HG, snoRNA distance to upstream and downstream exons, snoRNA
    distance to predicted branch_point, snoRNA structure stability (Minimal Free
    Energy or MFE), snoRNA terminal stem stability (MFE), snoRNA terminal
    stem length score, snoRNA conservation and hamming_distance_box (per box or
    global hamming distance). """

# Labels (abundance_cutoff and abundance_cutoff_2); feature abundance_cutoff_host
tpm_df_labels = pd.read_csv(snakemake.input.abundance_cutoff, sep='\t')
tpm_df_labels = tpm_df_labels[['gene_id', 'gene_name', 'abundance_cutoff', 'abundance_cutoff_2', 'abundance_cutoff_host']]
tpm_df_labels.columns = ['gene_id_sno', 'gene_name', 'abundance_cutoff', 'abundance_cutoff_2', 'abundance_cutoff_host']

# Other features
sno_length = pd.read_csv(snakemake.input.sno_length, sep='\t',
                names=['gene_id_sno', 'sno_length'])
snodb_nmd_di_promoters = pd.read_csv(snakemake.input.snodb_nmd_di_promoters, sep='\t')
snodb_nmd_di_promoters = snodb_nmd_di_promoters[['gene_id_sno', 'sno_type', 'sno_target',
                            'host_biotype2', 'NMD_susceptibility', 'di_promoter', 'host_function']]
sno_mfe = pd.read_csv(snakemake.input.sno_structure_mfe, sep='\t',
            names=['gene_id_sno', 'sno_mfe'])
terminal_stem_mfe = pd.read_csv(snakemake.input.terminal_stem_mfe, sep='\t',
                    names=['gene_id_sno', 'terminal_stem_mfe'])
terminal_stem_length_score = pd.read_csv(snakemake.input.terminal_stem_length_score,
                            sep='\t')
#conservation = pd.read_csv(snakemake.input.sno_conservation, sep='\t')
#conservation.columns = ['gene_id_sno', 'conservation_score']

hamming = pd.read_csv(snakemake.input.hamming_distance_box, sep='\t')
#hamming = hamming[['gene_id', 'C_hamming', 'D_hamming', 'C_prime_hamming', 'D_prime_hamming', 'H_hamming', 'ACA_hamming']]  # get hamming distance per box only (not combined)
hamming = hamming[['gene_id', 'combined_box_hamming']]  # get combined hamming distance for all boxes in a snoRNA
hamming = hamming.rename(columns={'gene_id': 'gene_id_sno'})

# Get sno location within intron (distances to branchpoint and to up/downstream exons)
location_bp = pd.read_csv(snakemake.input.location_and_branchpoint, sep='\t')
location_bp = location_bp[['gene_id_sno', 'intron_number', 'intron_length', 'exon_number_per_hg',
                        'distance_upstream_exon', 'distance_downstream_exon', 'dist_to_bp']]
# Change 'intron_number' col for 'intron_rank_5prime' (i.e. a better name) and compute actual total number of introns
# Compute also intron rank but counting from the 3' ('intron_rank_3prime') and relative_intron_rank (intron_rank_5prime / intron_number)
location_bp = location_bp.rename(columns={"intron_number": "intron_rank_5prime"})
location_bp['total_intron_number'] = location_bp['exon_number_per_hg'] - 1
location_bp['intron_rank_3prime'] = location_bp['exon_number_per_hg'] - location_bp['intron_rank_5prime']
location_bp.loc[location_bp['intron_rank_5prime'] == 0, 'intron_rank_3prime'] = 0  # patch for snoRNAs encoded (completely or with an overlap) within an exon of a HG
location_bp['relative_intron_rank'] = location_bp['intron_rank_3prime'] / location_bp['total_intron_number']
location_bp = location_bp.replace(np.inf, 0)  # replace relative_intron_rank to 0 when total_intron_number is equal to  0
location_bp = location_bp.drop(columns=['exon_number_per_hg'])


# Merge iteratively all of these dataframes
#df_list = [tpm_df_labels, sno_length, conservation, hamming, snodb_nmd_di_promoters,
 #           location_bp, sno_mfe, terminal_stem_mfe, terminal_stem_length_score]
df_list = [tpm_df_labels, sno_length, hamming, snodb_nmd_di_promoters,
            location_bp, sno_mfe, terminal_stem_mfe, terminal_stem_length_score]

df_label = df_list[0]
temp = [df_label]
for i, df in enumerate(df_list[1:]):
    if i == 0:
        df_temp = temp[0].merge(df, how='left', on='gene_id_sno')
        temp.append(df_temp)
    else:
        df_temp = temp[i].merge(df, how='left', on='gene_id_sno')
        temp.append(df_temp)

final_df = temp[-1]  # get the last df in temp, i.e. the final merged df

final_df.to_csv(snakemake.output.feature_df, index=False, sep='\t')
