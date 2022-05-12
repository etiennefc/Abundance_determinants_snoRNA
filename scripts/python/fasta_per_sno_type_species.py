#!/usr/bin/python3
import pandas as pd

sno_fasta = snakemake.input.sno_fasta
sno_info = pd.read_csv(snakemake.input.sno_info, sep='\t')
cd_output, haca_output = snakemake.output.cd_fasta, snakemake.output.haca_fasta

# Select either C/D or H/ACA box snoRNAs in fasta of all snoRNA sequences
cd_ids = sno_info[sno_info['snoRNA_type'] == 'C/D']['gene_id'].to_list()
haca_ids = sno_info[sno_info['snoRNA_type'] == 'H/ACA']['gene_id'].to_list()
cd_dict, haca_dict = {}, {}
with open(sno_fasta, 'r') as f:
    sno_id = ''
    for line in f:
        if line.startswith('>'):
                id = line.lstrip('>').rstrip('\n')
                sno_id = id
        else:
            seq = line.rstrip('\n')
            if sno_id in cd_ids:
                cd_dict[sno_id] = seq
            elif sno_id in haca_ids:
                haca_dict[sno_id] = seq

# Create fasta of C/D snoRNAs
with open(cd_output, 'w') as f:
    for sno_id, sequence in cd_dict.items():
        f.write(f'>{sno_id}\n')
        f.write(f'{sequence}\n')

# Create fasta of H/ACA snoRNAs
with open(haca_output, 'w') as f:
    for sno_id, sequence in haca_dict.items():
        f.write(f'>{sno_id}\n')
        f.write(f'{sequence}\n')



