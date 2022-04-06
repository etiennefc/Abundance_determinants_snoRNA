#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool

""" Get sno sequences from snoRNA bed file and genome fasta."""

bed = BedTool(snakemake.input.sno_bed)
genome = snakemake.input.genome
output = snakemake.output.sno_fasta

fasta = bed.sequence(fi=genome, nameOnly=True, rna=True)
with open(fasta.seqfn, 'r') as fasta_file, open(output, 'w') as output_file:
    for line in fasta_file:
        output_file.write(line)
