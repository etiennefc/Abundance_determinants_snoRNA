#!/usr/bin/python3
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
import subprocess as sp

input_file_stem = snakemake.input.snora77b_terminal_stem
input_file_all_sno = snakemake.input.all_snorna_structure
snora77b_dot_bracket = snakemake.output.snora77b_dot_bracket

# Remove temporarily the MFE (i.e a negative number between parentheses) in the snoRNA/terminal stem stability file
sp.call(f"sed -E 's/\(.[0-9]*.[0-9]*\)//g' {input_file_stem} > temp_stem.fa", shell=True)

# Extract the sequence and dot bracket of SNORA77B only from the fasta of all sno dot brackets
sp.call(f"grep -A 2 'ENSG00000264346' {input_file_all_sno} > {snora77b_dot_bracket}", shell=True)
sp.call(f"sed -E 's/\(.[0-9]*.[0-9]*\)//g' {snora77b_dot_bracket} > temp_sno.fa", shell=True)

# Create forgi graph of SNORA77B with its terminal stem
plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
cg = forgi.load_rna("temp_stem.fa", allow_many=False)
fvm.plot_rna(cg, text_kwargs={"fontweight":"black"}, lighten=0.7,
             backbone_kwargs={"linewidth":2}, ax=ax)
plt.savefig(snakemake.output.snora77b_terminal_stem_figure)

# Create forgi graph of SNORA77B alone
plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
cg = forgi.load_rna("temp_sno.fa", allow_many=False)
fvm.plot_rna(cg, text_kwargs={"fontweight":"black"}, lighten=0.7,
             backbone_kwargs={"linewidth":2}, ax=ax)
plt.savefig(snakemake.output.snora77b_figure)

# Remove temp files
sp.call('rm temp_s*', shell=True)
