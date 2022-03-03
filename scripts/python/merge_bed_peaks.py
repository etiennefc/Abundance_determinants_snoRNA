#!/usr/bin/python3
from pybedtools import BedTool

unmerged_beds = snakemake.input.input_beds
output_beds = snakemake.output.merge_beds

# Merge overlapping peaks for each bed file by summing the score value of overlapping peaks
for i, path in enumerate(unmerged_beds):
    bed = BedTool(path)
    merged_bed = bed.merge(s=True, c="5,6", o="sum,distinct").saveas(output_beds[i])
