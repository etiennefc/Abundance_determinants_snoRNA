import os

include: "data_processing.smk"
include: "downloads.smk"

rule flank_extend_snoRNA:
    """ Get the upstream and downstream flanking regions (15 nt) of each snoRNA
        using pybedtools flank. Then extend these flanking regions inside the
        snoRNA sequence using pybedtools slop. Extend for 5 nt from the 5' and
        3' of C/D box snoRNAs; extend 5 nt from the 5' and 3 nt from the 3' of
        H/ACA box snoRNAs."""
    input:
        all_sno_bed = rules.gtf_to_bed.output.all_sno_bed,
        snodb_info = rules.format_snodb.output.snodb_formatted
    output:
        flanking_cd_left = config['path']['flanking_regions_cd_left'],
        flanking_cd_right = config['path']['flanking_regions_cd_right'],
        flanking_haca_left = config['path']['flanking_regions_haca_left'],
        flanking_haca_right = config['path']['flanking_regions_haca_right']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/flank_extend_snoRNA.py"


rule add_chr_to_genome:
    """ Add the prefix 'chr' before chromosome number in the downloaded genome
        fasta file."""
    input:
        genome = rules.ensembl_genome.output.genome
    output:
        genome_chr = config['path']['genome_chr']
    shell:
        "sed 's/>/>chr/g' {input.genome} > {output.genome_chr}"


rule get_fasta:
    """ Get the sequence of all the extended flanking regions into a fasta file
        using pybedtools 'sequence'."""
    input:
        genome_fasta = rules.add_chr_to_genome.output.genome_chr,
        flanking_cd_left = rules.flank_extend_snoRNA.output[0],
        flanking_cd_right = rules.flank_extend_snoRNA.output[1],
        flanking_haca_left = rules.flank_extend_snoRNA.output[2],
        flanking_haca_right = rules.flank_extend_snoRNA.output[3]
    output:
        sequences = config['path']['flanking_regions_sno_fasta']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_fasta_terminal_stem.py"


rule rna_cofold:
    """ Get the structure stability and number of base paired in the potential
        terminal stem of snoRNAs."""
    input:
        fasta = rules.get_fasta.output.sequences
    output:
        mfe_stem = config['path']['terminal_stem_mfe_fasta']
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "RNAcofold < {input.fasta} > {output.mfe_stem} && "
        "mv *.ps data/terminal_stem/"


rule fasta_to_tsv_terminal_stem_mfe:
    """ Convert the fasta output of RNAcofold into a tsv table with a snoRNA id
        as the first column and the MFE of the terminal stem as the second column."""
    input:
        mfe = rules.rna_cofold.output.mfe_stem
    params:
        temp_mfe = "temporary_mfe.tsv",
        temp_id = "temporary_ids.tsv"
    output:
        mfe_stem_final = config['path']['terminal_stem_mfe_tsv']
    shell:
        """grep -E ">" {input.mfe} | sed 's/>//g' > {params.temp_id} && """
        """grep -oE "\-*[0-9]+\.[0-9]*" {input.mfe} > {params.temp_mfe} && """
        """paste {params.temp_id} {params.temp_mfe} > {output.mfe_stem_final} && """
        """rm temporary_*"""


rule get_terminal_stem_length:
    """ Get the length of the potential terminal stems associated to each snoRNA. """
    input:
        rna_cofold = rules.rna_cofold.output.mfe_stem
    output:
        length_stem = config['path']['terminal_stem_length']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/terminal_stem_length.py"
