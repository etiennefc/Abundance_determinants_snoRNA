#!/usr/bin/python3
import pandas as pd
import re
""" Get the terminal stem length score of snoRNAs, if they have a realistic
    terminal stem. The score is equal to the intermolecular paired nt between
    the left and right flanking regions of snoRNAs minus the number of nt gaps
    within the stem."""

dot_bracket_mfe_fasta = snakemake.input.rna_cofold

sno_id = ['']
terminal_stem = {}
with open(dot_bracket_mfe_fasta, 'r') as file:
    for line in file:
        if line.startswith('>'):  # sno_id lines
            sno_id_clean = str(line)
            sno_id[0] = sno_id_clean[1:].replace('\n', '')
        elif '(' in line:  #dot_bracket lines
            dot_bracket = str(line)
            dot_bracket = dot_bracket[0:20].replace('\n', '')  # we select only the left flanking region (20 nt)
            paired_base = dot_bracket.count('(')
            intramolecular_paired = dot_bracket.count(')')  # these ')' are intramolecular paired nt
                                                            # (i.e. nt pair within left sequence only)
            # This finds all overlapping (and non-overlapping) gaps of 1 to 19 nt inside the left flanking region
            gaps = re.findall(r'(?=(\(\.{1,19}\())', dot_bracket)
            number_gaps = ''.join(gaps)  # join all gaps together in one string
            number_gaps = len(re.findall('\.', number_gaps))  # count the number of nt gaps in sequence
            #print(sno_id[0])
            #print(dot_bracket)
            #print(gaps)
            #print(number_gaps)
            stem_length_score = paired_base - intramolecular_paired - number_gaps
            if stem_length_score < 0:
                stem_length_score = 0
            terminal_stem[sno_id[0]] = stem_length_score

df = pd.DataFrame.from_dict(terminal_stem, orient='index').reset_index()
df.columns = ['gene_id_sno', 'terminal_stem_length_score']

df.to_csv(snakemake.output.length_stem, index=False, sep='\t')
