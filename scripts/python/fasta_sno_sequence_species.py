#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

""" Get sno sequences from snoRNA bed file and genome fasta."""

bed = BedTool(snakemake.input.sno_bed)
genome = snakemake.input.genome
if 'Mus_musculus' in genome:
    temp_output = 'mouse_temp_output.fa'
else:
    species = snakemake.wildcards.species
    temp_output = f'{species}_temp_output.fa'
output = snakemake.output.sno_fasta

fasta = bed.sequence(fi=genome, nameOnly=True, s=True)
with open(fasta.seqfn, 'r') as fasta_file, open(temp_output, 'w') as output_file:
    for line in fasta_file:
        if '>' not in line:
            line = line.replace('T', 'U')
        output_file.write(line)

# Remove strand info from fasta ids
sp.call(f"sed -i 's/(+)//g; s/(-)//g' {temp_output} && mv {temp_output} {output}", shell=True)
