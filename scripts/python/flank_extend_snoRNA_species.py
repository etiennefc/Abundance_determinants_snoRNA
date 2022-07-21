#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

""" Get the upstream and downstream flanking regions (15 nt) of each snoRNA
    using pybedtools flank. Then extend these flanking regions inside the
    snoRNA sequence using pybedtools slop. Extend for 5 nt from the 5' and
    3' of C/D box snoRNAs; extend 5 nt from the 5' and 3 nt from the 3' of
    H/ACA box snoRNAs."""

species = snakemake.wildcards.species
all_sno_bed = BedTool(snakemake.input.all_sno_bed)
sno_info_df = pd.read_csv(snakemake.input.sno_info, sep='\t')
chr_size_file = snakemake.input.genome_chr_size

# Get snoRNA type of all snoRNAs from RNAcentral reference
sno_type_dict = sno_info_df.set_index('gene_id')['snoRNA_type'].to_dict()

# Filter sno_bed to only keep snoRNAs where we could find a snoRNA type in RNAcentral
all_sno_bed = BedTool(line for line in all_sno_bed if line[3] in sno_type_dict.keys())

# Get 15 nt flanking regions upstream and downstream of C/D and H/ACA snoRNAs
# The .saveas() is needed to create a temporary version of the object since it's used afterwards (otherwise, it does not work)
flanking = all_sno_bed.flank(g=chr_size_file, b=15) #this is a bedtools object


cd_flank = BedTool(line for line in flanking if sno_type_dict[line[3]] == 'C/D').saveas()  # where line[3] corresponds to the gene_id
haca_flank = BedTool(line for line in flanking if sno_type_dict[line[3]] == 'H/ACA').saveas()  # where line[3] corresponds to the gene_id

# Separate the bed objects into the left or the right flanking regions for both snoRNA type
cd_flank_left = BedTool(line for i, line in enumerate(cd_flank) if i % 2 == 0).saveas()
cd_flank_right = BedTool(line for i, line in enumerate(cd_flank) if i % 2 != 0).saveas()

haca_flank_left = BedTool(line for i, line in enumerate(haca_flank) if i % 2 == 0).saveas()
haca_flank_right = BedTool(line for i, line in enumerate(haca_flank) if i % 2 != 0).saveas()

# For H/ACA snoRNAs, split by strand, since we don't extend (slop) the same number of nt (5 vs 3 nt respectively for the 5' and 3' end)
# We don't need to that for C/D since we extend 5 nt from the 5' and 3'end, so it does not affect the terminal stem sequences
haca_flank_left_plus = BedTool(line for line in haca_flank_left if line[5] == '+').saveas()  # where line[5] corresponds to the strand
haca_flank_left_minus = BedTool(line for line in haca_flank_left if line[5] == '-').saveas()  # where line[5] corresponds to the strand
haca_flank_right_plus = BedTool(line for line in haca_flank_right if line[5] == '+').saveas()  # where line[5] corresponds to the strand
haca_flank_right_minus = BedTool(line for line in haca_flank_right if line[5] == '-').saveas()  # where line[5] corresponds to the strand

# For C/D snoRNAs, extend the flanking region inside the snoRNA for 5 nt from the 5' and 3' of the snoRNA
cd_extend_left = cd_flank_left.slop(r=5, l=0, g=chr_size_file).saveas(snakemake.output.flanking_cd_left)
cd_extend_right = cd_flank_right.slop(l=5, r=0, g=chr_size_file).saveas(snakemake.output.flanking_cd_right)

# For H/ACA snoRNAs, extend the flanking region inside the snoRNA for 5 nt from the 5' and 3 nt from the 3' of the snoRNA
haca_extend_left_plus = haca_flank_left_plus.slop(r=5, l=0, g=chr_size_file, s=True).saveas(f'temp_haca_extend_left_plus_{species}.bed')
haca_extend_right_plus = haca_flank_right_plus.slop(l=3, r=0, g=chr_size_file, s=True).saveas(f'temp_haca_extend_right_plus_{species}.bed')

haca_extend_left_minus = haca_flank_left_minus.slop(l=3, r=0, g=chr_size_file, s=True).saveas(f'temp_haca_extend_left_minus_{species}.bed')
haca_extend_right_minus = haca_flank_right_minus.slop(r=5, l=0, g=chr_size_file, s=True).saveas(f'temp_haca_extend_right_minus_{species}.bed')

sp.call(f'cat temp_haca_extend_left_plus_{species}.bed temp_haca_extend_left_minus_{species}.bed > '+snakemake.output.flanking_haca_left, shell=True)
sp.call(f'cat temp_haca_extend_right_plus_{species}.bed temp_haca_extend_right_minus_{species}.bed > '+snakemake.output.flanking_haca_right, shell=True)
sp.call(f'rm temp_haca_extend*_{species}.bed', shell=True)

