#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

""" Determine the extended overlap between a bed file of all snoRNAs and a bed
    of the binding of AQR (eCLIP data)."""

col = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature',
        'dot2', 'gene_info']
sno_bed = pd.read_csv(snakemake.input.sno_bed, sep='\t', names=col)  # generated with gtf_to_bed
HG_bed = BedTool(snakemake.input.HG_bed)
aqr_hepg2, aqr_k562 = snakemake.input.aqr_HepG2_bed, snakemake.input.aqr_K562_bed
df = pd.read_csv(snakemake.input.df, sep='\t')


# First, merge the two AQR eCLIP bed files (remove redundancy in peaks)
# All peaks have at least a pVal < 0.01 and max 1 nt between them
sp.call(f'cat {aqr_hepg2} {aqr_k562} | sort -k1,1 -k2,2n > temp_aqr.bed', shell=True)
aqr_temp_bed = BedTool('temp_aqr.bed')
aqr_bed = aqr_temp_bed.merge(s=True, d=1, c=[5,7,6,8], o='distinct')


# Generate bed of the intron of all intronic snoRNAs
df = df[df['abundance_cutoff_host'] != 'intergenic']
sno_bed = sno_bed[sno_bed['gene_id'].isin(list(df.gene_id_sno))]
d_upstream = dict(zip(df.gene_id_sno, df.distance_upstream_exon))
d_downstream = dict(zip(df.gene_id_sno, df.distance_downstream_exon))
intron_bed = sno_bed.copy()

# Extend sno_bed to upstream and downstream exon (thereby the snoRNA whole intron) depending on the strand
intron_bed.loc[intron_bed['strand'] == "+", 'start_intron'] = intron_bed['start'] - intron_bed['gene_id'].map(d_upstream)
intron_bed.loc[intron_bed['strand'] == "+", 'end_intron'] = intron_bed['end'] + intron_bed['gene_id'].map(d_downstream)
intron_bed.loc[intron_bed['strand'] == "-", 'start_intron'] = intron_bed['start'] - intron_bed['gene_id'].map(d_downstream)
intron_bed.loc[intron_bed['strand'] == "-", 'end_intron'] = intron_bed['end'] + intron_bed['gene_id'].map(d_upstream)
intron_bed = intron_bed[['chr', 'start_intron', 'end_intron', 'gene_id', 'dot', 'strand', 'source', 'feature',
        'dot2', 'gene_info']]
intron_bed['start_intron'] = intron_bed['start_intron'].astype('int')
intron_bed['end_intron'] = intron_bed['end_intron'].astype('int')

intron_bed.to_csv(snakemake.output.intron_bed, sep='\t', header=False, index=False)
intron_bed = BedTool(snakemake.output.intron_bed)

# Intersect the aqr peaks with snoRNA intron bed
intersection = intron_bed.intersect(aqr_bed, wa=True, s=True, wb=True, sorted=True).saveas(snakemake.output.overlap_sno_AQR)

sp.call('rm temp_aqr.bed', shell=True)

