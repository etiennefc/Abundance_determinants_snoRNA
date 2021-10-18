#!/usr/bin/python3
import pandas as pd

""" Generate a merged dataframe of all snoRNA features and labels that will be
    used by the predictor. The labels are 'expressed' or 'not_expressed' using
    the column 'abundance_cutoff_2'. The features are: snoRNA length, type,
    target, host gene (HG) biotype and function, HG NMD susceptibility, HG
    promoter type, intron number and length in which the snoRNA is encoded,
    snoRNA distance to upstream and downstream exons, snoRNA distance to
    predicted branch_point, snoRNA structure stability (Minimal Free Energy or
    MFE), snoRNA terminal stem stability (MFE) and snoRNA terminal stem length
    score. """

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
location_bp = pd.read_csv(snakemake.input.location_and_branchpoint, sep='\t')
location_bp = location_bp[['gene_id_sno', 'intron_number', 'intron_length',
                        'distance_upstream_exon', 'distance_downstream_exon', 'dist_to_bp']]
sno_mfe = pd.read_csv(snakemake.input.sno_structure_mfe, sep='\t',
            names=['gene_id_sno', 'sno_mfe'])
terminal_stem_mfe = pd.read_csv(snakemake.input.terminal_stem_mfe, sep='\t',
                    names=['gene_id_sno', 'terminal_stem_mfe'])
terminal_stem_length_score = pd.read_csv(snakemake.input.terminal_stem_length_score,
                            sep='\t')
conservation = pd.read_csv(snakemake.input.sno_conservation, sep='\t')
conservation.columns = ['gene_id_sno', 'conservation_score']

# Merge iteratively all of these dataframes
df_list  = [tpm_df_labels, sno_length, conservation, snodb_nmd_di_promoters,
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
