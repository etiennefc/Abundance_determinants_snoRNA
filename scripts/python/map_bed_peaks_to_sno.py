#!/usr/bin/python3
from pybedtools import BedTool

merged_beds = snakemake.input.merged_beds
snoRNA_bed = BedTool(snakemake.input.snoRNA_bed)
output_beds = snakemake.output.mapped_snoRNA_bed

# Map the enrichment value of overlapping peaks to each snoRNA (return the sum of all peaks
#overlapping a snoRNA per RBP).
for i, path in enumerate(merged_beds):
    bed = BedTool(path)
    mapped_bed = snoRNA_bed.map(b=bed, s=True, c="4", o="sum").saveas(output_beds[i])
