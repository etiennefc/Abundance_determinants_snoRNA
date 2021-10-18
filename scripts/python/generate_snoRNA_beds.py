#!/usr/bin/python3
import pandas as pd
import subprocess as sp
""" From a bed file containing all snoRNAs (generated via gtf_to_bed.sh) and a csv file containing snoRNA info (i.e.
    their HG), generate snoRNA bed files (intronic, intergenic, intronic without SNHG14 snoRNAs, only SNHG14 snoRNAs)"""

all_sno_bed = pd.read_csv(snakemake.input.all_sno_bed, sep='\t',
                          names=['chr', 'start', 'end', 'gene_id', 'dot',
                                'strand', 'source', 'feature', 'dot2',
                                 'gene_info'])
all_sno_bed['chr'] = all_sno_bed['chr'].astype(str)
all_sno_bed = all_sno_bed.sort_values(['chr', 'start', 'end'])
all_sno_bed = all_sno_bed.replace(to_replace='"cluster_([0-9]{0,4})', value=r'"cluster_\1"; gene_version "1', regex=True)  # Add gene_version for blockbuster detected snoRNAs


sno_info = pd.read_csv(snakemake.input.sno_info_df)
intergenic_sno_id = list(sno_info[sno_info['host_id'].isna()].gene_id)  # snoRNAs with NaN host_id are intergenic
snhg14_sno_id = list(sno_info[sno_info['host_name'] == 'SNHG14'].gene_id)

intronic_sno_bed = all_sno_bed[~all_sno_bed['gene_id'].isin(intergenic_sno_id)]
intergenic_sno_bed = all_sno_bed[all_sno_bed['gene_id'].isin(intergenic_sno_id)]

intronic_sno_bed_wo_snhg14 = intronic_sno_bed[~intronic_sno_bed['gene_id'].isin(snhg14_sno_id)]  # bed file containing all intronic snoRNAs except those encoded in SNHG14
snhg14_sno_bed = intronic_sno_bed[intronic_sno_bed['gene_id'].isin(snhg14_sno_id)]  # bed file containing SNHG14 snoRNAs

intronic_sno_bed.to_csv(snakemake.output.intronic_sno_bed, index=False, header=False, sep='\t')
intergenic_sno_bed.to_csv(snakemake.output.intergenic_sno_bed, index=False, header=False, sep='\t')
intronic_sno_bed_wo_snhg14.to_csv(snakemake.output.intronic_sno_bed_wo_snhg14, index=False, header=False, sep='\t')
snhg14_sno_bed.to_csv(snakemake.output.snhg14_sno_bed, index=False, header=False, sep='\t')
