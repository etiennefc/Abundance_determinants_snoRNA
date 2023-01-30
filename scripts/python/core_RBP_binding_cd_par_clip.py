#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

""" Determine the overlap between a bed file of all C/D snoRNAs and a bed
    of the binding of different C/D core proteins (PAR-CLIP data)."""

col = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature',
        'dot2', 'gene_info']
sno_bed = pd.read_csv(snakemake.input.sno_bed, sep='\t', names=col)  # generated with gtf_to_bed
nop58_a_bed = [path for path in snakemake.input.RBP_par_clip_beds if 'NOP58_repA' in path][0]
nop58_b_bed = [path for path in snakemake.input.RBP_par_clip_beds if 'NOP58_repB' in path][0]
nop56_bed = [path for path in snakemake.input.RBP_par_clip_beds if 'NOP56' in path][0]
fbl_bed = [path for path in snakemake.input.RBP_par_clip_beds if 'FBL_merge' in path][0]
fbl_mnase_bed = [path for path in snakemake.input.RBP_par_clip_beds if 'FBL_mnase' in path][0]
df = pd.read_csv(snakemake.input.df, sep='\t')

# Keep only C/D snoRNAs in sno_bed 
cd_bed = sno_bed[sno_bed['gene_id'].isin(list(df[df['sno_type'] == 'C/D'].gene_id_sno))]
cd_bed.to_csv('cd_temp.bed', index=False, header=False, sep='\t')
cd_bed = BedTool('cd_temp.bed')

# Merge the RBP PAR-CLIP bed file (remove redundancy in peaks)
# Merged peaks have at most 1 nt between them. Concat the 2 NOP58 replicates together (same for FBL)
sp.call(f'cat {nop58_a_bed} {nop58_b_bed} | sort -k1,1 -k2,2n > nop58_temp_merge.bed', shell=True)
nop58_parclip_temp_bed = BedTool('nop58_temp_merge.bed')
nop58_parclip_bed = nop58_parclip_temp_bed.merge(s=True, d=1, c=[6,5,6], o='distinct,sum,distinct').saveas()

sp.call(f'cat {fbl_bed} {fbl_mnase_bed} | sort -k1,1 -k2,2n > fbl_temp_merge.bed', shell=True)
fbl_parclip_temp_bed = BedTool('fbl_temp_merge.bed')
fbl_parclip_bed = fbl_parclip_temp_bed.merge(s=True, d=1, c=[6,5,6], o='distinct,sum,distinct').saveas()

nop56_parclip_temp_bed = BedTool(nop56_bed)
nop56_parclip_bed = nop56_parclip_temp_bed.merge(s=True, d=1, c=[6,5,6], o='distinct,sum,distinct').saveas()

# Intersect the RBP peaks with C/D snoRNA bed (make sure that least 50% of the RBP peak is overlapped by a given snoRNA)
intersection = cd_bed.intersect(nop58_parclip_bed, wa=True, s=True, wb=True, sorted=True, F=0.5).saveas(snakemake.output.overlap_sno_NOP58_par_clip)
intersection2 = cd_bed.intersect(fbl_parclip_bed, wa=True, s=True, wb=True, sorted=True, F=0.5).saveas(snakemake.output.overlap_sno_FBL_par_clip)
intersection3 = cd_bed.intersect(nop56_parclip_bed, wa=True, s=True, wb=True, sorted=True, F=0.5).saveas(snakemake.output.overlap_sno_NOP56_par_clip)

sp.call('rm cd_temp.bed nop58_temp_merge.bed fbl_temp_merge.bed', shell=True)
