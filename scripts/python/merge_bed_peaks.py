#!/usr/bin/python3
from pybedtools import BedTool

unmerged_beds = snakemake.input.input_beds
unmerged_beds = [path for path in unmerged_beds if 'DKC1' not in path]
output_beds = snakemake.output.merge_beds
print(unmerged_beds)
print(output_beds)
# Merge overlapping peaks for each bed file by summing the score value of overlapping peaks (max 1 nt distance between peaks)
for i, path in enumerate(unmerged_beds):
    bed = BedTool(path)
    merged_bed = bed.merge(s=True, d=1, c="6,5,6", o="distinct,sum,distinct").saveas(output_beds[i])
