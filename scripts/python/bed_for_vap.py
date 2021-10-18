#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool

""" Generate specific bed files for expressed vs not_expressed snoRNAs (either
    intronic or intergenic) and of their HG (if the snoRNA is intronic)"""
cols = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature', 'dot2', 'characteristics']
intronic = pd.read_csv(snakemake.input.intronic_sno_bed, names=cols, sep='\t')
intergenic = pd.read_csv(snakemake.input.intergenic_sno_bed, names=cols, sep='\t')
hg_bed = pd.read_csv(snakemake.input.hg_bed, names=cols, sep='\t')
abundance_status_df = pd.read_csv(snakemake.input.abundance_status_df, sep='\t')
hg_df = pd.read_csv(snakemake.input.hg_df)

ab_status_dict = dict(zip(abundance_status_df['gene_id_sno'], abundance_status_df['abundance_cutoff_2']))

# Split intronic snoRNA bed file based on their abundance status
expressed_intronic_snoRNA = intronic.loc[intronic['gene_id'].isin([k for k,v in ab_status_dict.items() if v == 'expressed'])]
not_expressed_intronic_snoRNA = intronic.loc[intronic['gene_id'].isin([k for k,v in ab_status_dict.items() if v == 'not_expressed'])]

# Split intergenic snoRNA bed file based on their abundance status
expressed_intergenic_snoRNA = intergenic.loc[intergenic['gene_id'].isin([k for k,v in ab_status_dict.items() if v == 'expressed'])]
not_expressed_intergenic_snoRNA = intergenic.loc[intergenic['gene_id'].isin([k for k,v in ab_status_dict.items() if v == 'not_expressed'])]

# Create a bed file for HG of either expressed or not expressed snoRNAs
HG_expressed_sno_df = hg_df[hg_df['sno_id'].isin(list(expressed_intronic_snoRNA['gene_id']))]
HG_expressed_sno_bed = hg_bed[hg_bed['gene_id'].isin(list(HG_expressed_sno_df['host_id']))]

HG_not_expressed_sno_df = hg_df[hg_df['sno_id'].isin(list(not_expressed_intronic_snoRNA['gene_id']))]
HG_not_expressed_sno_bed = hg_bed[hg_bed['gene_id'].isin(list(HG_not_expressed_sno_df['host_id']))]

# Save all bed files
expressed_intronic_snoRNA.to_csv(snakemake.output.expressed_intronic_sno_bed, index=False, sep='\t', header=False)
not_expressed_intronic_snoRNA.to_csv(snakemake.output.not_expressed_intronic_sno_bed, index=False, sep='\t', header=False)
expressed_intergenic_snoRNA.to_csv(snakemake.output.expressed_intergenic_sno_bed, index=False, sep='\t', header=False)
not_expressed_intergenic_snoRNA.to_csv(snakemake.output.not_expressed_intergenic_sno_bed, index=False, sep='\t', header=False)
HG_expressed_sno_bed.to_csv(snakemake.output.HG_expressed_sno_bed, index=False, sep='\t', header=False)
HG_not_expressed_sno_bed.to_csv(snakemake.output.HG_not_expressed_sno_bed, index=False, sep='\t', header=False)
