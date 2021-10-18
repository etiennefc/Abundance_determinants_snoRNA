#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

""" Get the upstream and downstream flanking regions (15 nt) of each snoRNA
    using pybedtools flank. Then extend these flanking regions inside the
    snoRNA sequence using pybedtools slop. Extend for 5 nt from the 5' and
    3' of C/D box snoRNAs; extend 5 nt from the 5' and 3 nt from the 3' of
    H/ACA box snoRNAs."""


all_sno_bed = BedTool(snakemake.input.all_sno_bed)
sno_info_df = pd.read_csv(snakemake.input.snodb_info, sep='\t')

# Get snoRNA type of all snoRNAs from snoDB reference
sno_type_dict = sno_info_df.set_index('gene_id_sno')['sno_type'].to_dict()

# Get 15 nt flanking regions upstream and downstream of C/D and H/ACA snoRNAs
# The .saveas() is needed to create a temporary version of the object since it's used afeterwards (otherwise, it does not work)
flanking = all_sno_bed.flank(genome="hg38", b=15) #this is a bedtools object
cd_flank = BedTool(line for line in flanking if sno_type_dict[line[3]] == 'C/D').saveas()  # where line[3] corresponds to the gene_id
haca_flank = BedTool(line for line in flanking if sno_type_dict[line[3]] == 'H/ACA').saveas()  # where line[3] corresponds to the gene_id

# Separate the bed objects into the left or the right flanking regions for both snoRNA type
cd_flank_left = BedTool(line for i, line in enumerate(cd_flank) if i % 2 == 0).saveas()
cd_flank_right = BedTool(line for i, line in enumerate(cd_flank) if i % 2 != 0).saveas()

haca_flank_left = BedTool(line for i, line in enumerate(haca_flank) if i % 2 == 0).saveas()
haca_flank_right = BedTool(line for i, line in enumerate(haca_flank) if i % 2 != 0).saveas()


# For C/D snoRNAs, extend the flanking region inside the snoRNA for 5 nt from the 5' and 3' of the snoRNA
cd_extend_left = cd_flank_left.slop(r=5, l=0, genome="hg38").saveas(snakemake.output.flanking_cd_left)
cd_extend_right = cd_flank_right.slop(l=5, r=0, genome="hg38").saveas(snakemake.output.flanking_cd_right)

# For H/ACA snoRNAs, extend the flanking region inside the snoRNA for 5 nt from the 5' and 3 nt from the 3' of the snoRNA
haca_extend_left = haca_flank_left.slop(r=5, l=0, genome="hg38").saveas(snakemake.output.flanking_haca_left)
haca_extend_right = haca_flank_right.slop(l=3, r=0, genome="hg38").saveas(snakemake.output.flanking_haca_right)
