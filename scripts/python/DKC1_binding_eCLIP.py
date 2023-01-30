#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

""" Determine the overlap between a bed file of all H/ACA snoRNAs and a bed
    of the binding of DKC1 (eCLIP data)."""

col = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature',
        'dot2', 'gene_info']
sno_bed = pd.read_csv(snakemake.input.sno_bed, sep='\t', names=col)  # generated with gtf_to_bed
dkc1_hepg2, dkc1_hek293 = snakemake.input.dkc1_HepG2_eCLIP_bed, snakemake.input.dkc1_HEK293_par_clip_bed
df = pd.read_csv(snakemake.input.df, sep='\t')

# Keep only H/ACA snoRNAs in sno_bed (i.e. those that bind DKC1)
haca_bed = sno_bed[sno_bed['gene_id'].isin(list(df[df['sno_type'] == 'H/ACA'].gene_id_sno))]
haca_bed.to_csv('haca_temp.bed', index=False, header=False, sep='\t')
haca_bed = BedTool('haca_temp.bed')

# First, merge the DKC1 eCLIP bed file (remove redundancy in peaks)
# All peaks have at least a pVal < 0.01 and max 1 nt between them
dkc1_temp_bed = BedTool(dkc1_hepg2)
dkc1_bed = dkc1_temp_bed.merge(s=True, d=1, c=[5,7,6,8], o='distinct')

# Second, merge the DKC1 PAR-CLIP bed file (remove redundancy in peaks)
# All peaks have at least a pVal < 0.01 and max 1 nt between them
dkc1_parclip_temp_bed = BedTool(dkc1_hek293)
dkc1_parclip_bed = dkc1_parclip_temp_bed.merge(s=True, d=1, c=[6,5,6], o='distinct,sum,distinct').saveas()

# Intersect the DKC1 peaks with H/ACA snoRNA bed (make sure that least 50% of the DKC1 peak is overlapped by a given snoRNA)
intersection = haca_bed.intersect(dkc1_bed, wa=True, s=True, wb=True, sorted=True, F=0.5).saveas(snakemake.output.overlap_sno_DKC1_eCLIP)
intersection2 = haca_bed.intersect(dkc1_parclip_bed, wa=True, s=True, wb=True, sorted=True, F=0.5).saveas(snakemake.output.overlap_sno_DKC1_par_clip)
sp.call('rm haca_temp.bed', shell=True)
