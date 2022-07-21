import os

include: "data_processing.smk"

rule format_snoRNA_sequences:
    """ Format snoRNA DNA sequence into RNA and create a fasta out of it."""
    input:
        snodb = rules.format_snodb.output.snodb_formatted
    output:
        sno_sequences = config['path']['sno_sequences']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/clean_sno_sequences.py"

rule rna_fold:
    """ Get the structure stability of snoRNAs (their MFE or Minimal Free
        Energy) using RNAfold."""
    input:
        sequences = rules.format_snoRNA_sequences.output.sno_sequences
    params:
        temp_name = "sno.mfe"
    output:
        mfe = config['path']['structure_stability_fasta']
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "sed 's/>/>HUMAN_/g' {input.sequences} > HUMAN_rna_fold.tsv && "
        "RNAfold --infile=HUMAN_rna_fold.tsv --outfile={params.temp_name} && "
        "sed -i 's/HUMAN_//g' {params.temp_name} && "
        "mv {params.temp_name} {output.mfe} && "
        "mv HUMAN_*.ps data/structure/stability/ && "
        "rm -f data_structure_stability_mfe.tsv HUMAN_rna_fold.tsv"

rule fasta_to_tsv:
    """ Convert the fasta output of RNA fold into a tsv table with a snoRNA id
        as the first column and the MFE as the second column."""
    input:
        mfe = rules.rna_fold.output.mfe
    params:
        temp_mfe = "temporary_mfe.tsv",
        temp_id = "temporary_ids.tsv"
    output:
        mfe_final = config['path']['structure_stability_tsv']
    shell:
        """grep -E ">" {input.mfe} | sed 's/>//g' > {params.temp_id} && """
        """grep -oE "\-*[0-9]+\.[0-9]*" {input.mfe} > {params.temp_mfe} && """
        """paste {params.temp_id} {params.temp_mfe} > {output.mfe_final} && """
        """rm temporary_*"""

rule RNAalifold_snora77b_stem:
    """ Get consensus sequence of SNORA77B+terminal stem aligned sequences in
        all vertebrates using as input the stockholm (stk) mulitple alignment
        format."""
    input:
        stk = "30primates_mod.stk"
    output:
        consensus_sequence = "CS.stk"
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "RNAalifold -f S {input.stk} > {output.consensus_sequence}"
