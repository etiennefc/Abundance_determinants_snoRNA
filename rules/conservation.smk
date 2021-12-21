import os

include: "downloads.smk"
include: "data_processing.smk"
include: "structure.smk"

rule bw_to_bg:
    """ Convert phastCons bigwig file into bedgraph file and sort that bedgraph
        in place by chromosome and start."""
    input:
        phastcons_bigwig = rules.phastcons_download.output.phastcons
    output:
        phastcons_bedgraph = config['path']['phastcons_bg']
    conda:
        "../envs/conservation.yaml"
    shell:
        """bigWigToBedGraph {input.phastcons_bigwig} {output.phastcons_bedgraph} """

rule sno_conservation:
    """ For each snoRNA, get the conservation score (phastCons) across a 100
        vertebrates for every nt composing a snoRNA (0 being not conserved and
        1 being conserved in all vertebrates). Then take the average score
        across all nt to generate a conservation score for each snoRNA."""
    input:
        phastcons_bg = rules.bw_to_bg.output.phastcons_bedgraph,
        sno_bed = rules.gtf_to_bed.output.all_sno_bed
    output:
        intersection_sno_conservation = config['path']['intersection_sno_conservation'],
        sno_conservation = config['path']['sno_conservation']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/sno_conservation.py"

rule fasta_per_sno_type:
    """ Create a fasta per snoRNA type out of foasta containing all snoRNAs."""
    input:
        sno_fasta = rules.format_snoRNA_sequences.output.sno_sequences,
        snodb = rules.format_snodb.output.snodb_formatted
    output:
        cd_fasta = config['path']['all_cd_fa'],
        haca_fasta = config['path']['all_haca_fa']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/fasta_per_sno_type.py"


rule c_d_box_location_all:
    """ Find C, D, C' and D' motif and location (if they exist) in ALL C/D box
        snoRNAs. The output of this rule will be used to create the hamming
        distance feature for all boxes of C/D snoRNAs in the initial feature_df."""
    input:
        cd_fasta = rules.fasta_per_sno_type.output.cd_fasta
    output:
        c_d_box_location = config['path']['c_d_and_prime_box_location']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/cd_box_location_all.py"


rule h_aca_box_location_all:
    """ Find H and ACA motif and location (if they exist) in ALL H/ACA box
        snoRNAs. The output of this rule will be used to create the hamming
        distance feature for all boxes of H/ACA snoRNAs in the initial feature_df."""
    input:
        dot_bracket = rules.rna_fold.output.mfe,
        haca_fasta = rules.fasta_per_sno_type.output.haca_fasta
    output:
        h_aca_box_location = config['path']['h_aca_box_location']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/haca_box_location_all.py"

rule hamming_distance_box_all:
    """ Get the hamming distance per box and also a combined hamming distance
        per snoRNA boxes for all snoRNAs as input features of the models."""
    input:
        c_d_box_location = rules.c_d_box_location_all.output.c_d_box_location,
        h_aca_box_location = rules.h_aca_box_location_all.output.h_aca_box_location
    output:
        hamming_distance_box_df = config['path']['hamming_distance_box_df']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/hamming_distance_box_all.py"


rule fasta_sequence_abundance_status:
    """ Generate a fasta of snoRNA sequence per sno_type and per abundance status
        (that will be useful for consensus_sequence_abundance_status)."""
    input:
        sno_fasta = config['path']['sno_sequences'],
        all_features_labels_df = config['path']['feature_df']
    output:
        expressed_cd = config['path']['expressed_cd_fa'],
        expressed_haca = config['path']['expressed_haca_fa'],
        not_expressed_cd = config['path']['not_expressed_cd_fa'],
        not_expressed_haca = config['path']['not_expressed_haca_fa']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/fasta_sequence_abundance_status.py"

rule fasta_sequence_abundance_status_length:
    """ Split the snoRNA fasta for C/D box snoRNAs into subgroups according to
        their length to see if they better align (below (small) or above (long)
        200 nt)."""
    input:
        sno_fasta = config['path']['sno_sequences'],
        all_features_labels_df = config['path']['feature_df']
    output:
        small_expressed_cd = config['path']['small_expressed_cd_fa'],
        long_expressed_cd = config['path']['long_expressed_cd_fa'],
        small_not_expressed_cd = config['path']['small_not_expressed_cd_fa'],
        long_not_expressed_cd = config['path']['long_not_expressed_cd_fa']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/fasta_sequence_abundance_status_length.py"

rule c_d_box_location:
    """ Find C, D, C' and D' motif and location (if they exist) in expressed vs
        not expressed C/D snoRNAs."""
    input:
        expressed_cd = rules.fasta_sequence_abundance_status.output.expressed_cd,
        not_expressed_cd = rules.fasta_sequence_abundance_status.output.not_expressed_cd
    output:
        c_d_box_location_expressed = config['path']['c_d_and_prime_box_location_expressed'],
        c_d_box_location_not_expressed = config['path']['c_d_and_prime_box_location_not_expressed']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/cd_box_location.py"

rule h_aca_location:
    """ Find H and ACA motif (if they exist in) in expressed and not expressed
        H/ACA snoRNAs."""
    input:
        dot_bracket = rules.rna_fold.output.mfe,
        expressed_haca = rules.fasta_sequence_abundance_status.output.expressed_haca,
        not_expressed_haca = rules.fasta_sequence_abundance_status.output.not_expressed_haca
    output:
        h_aca_box_location_expressed = config['path']['h_aca_box_location_expressed'],
        h_aca_box_location_not_expressed = config['path']['h_aca_box_location_not_expressed']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/haca_box_location.py"

rule flanking_nt_to_motifs:
    """ Get 3 nt flanking both up- and downstream of all motifs and create a
        fasta of motif plus flanking nucleotides per abundance status and sno_type."""
    input:
        expressed_cd_fa = rules.fasta_sequence_abundance_status.output.expressed_cd,
        not_expressed_cd_fa = rules.fasta_sequence_abundance_status.output.not_expressed_cd,
        expressed_haca_fa = rules.fasta_sequence_abundance_status.output.expressed_haca,
        not_expressed_haca_fa = rules.fasta_sequence_abundance_status.output.not_expressed_haca,
        c_d_box_location_expressed = rules.c_d_box_location.output.c_d_box_location_expressed,
        c_d_box_location_not_expressed = rules.c_d_box_location.output.c_d_box_location_not_expressed,
        h_aca_box_location_expressed = rules.h_aca_location.output.h_aca_box_location_expressed,
        h_aca_box_location_not_expressed = rules.h_aca_location.output.h_aca_box_location_not_expressed
    output:
        c_expressed = config['path']['c_box_expressed_flanking_nt_fa'],
        c_not_expressed = config['path']['c_box_not_expressed_flanking_nt_fa'],
        d_expressed = config['path']['d_box_expressed_flanking_nt_fa'],
        d_not_expressed = config['path']['d_box_not_expressed_flanking_nt_fa'],
        c_prime_expressed = config['path']['c_prime_box_expressed_flanking_nt_fa'],
        c_prime_not_expressed = config['path']['c_prime_box_not_expressed_flanking_nt_fa'],
        d_prime_expressed = config['path']['d_prime_box_expressed_flanking_nt_fa'],
        d_prime_not_expressed = config['path']['d_prime_box_not_expressed_flanking_nt_fa'],
        h_expressed = config['path']['h_box_expressed_flanking_nt_fa'],
        h_not_expressed = config['path']['h_box_not_expressed_flanking_nt_fa'],
        aca_expressed = config['path']['aca_box_expressed_flanking_nt_fa'],
        aca_not_expressed = config['path']['aca_box_not_expressed_flanking_nt_fa']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/flanking_nt_to_motifs.py"

rule convert_motif_table_to_fa:
    """ Convert table output of c_d_box_location and h_aca_location to fasta per
        motif (ex: fasta of all D box sequences per expressed C/D snoRNAs). The
        output fastas will be used to construct logo of motifs per abundance status"""
    input:
        c_d_box_location_expressed = rules.c_d_box_location.output.c_d_box_location_expressed,
        c_d_box_location_not_expressed = rules.c_d_box_location.output.c_d_box_location_not_expressed,
        h_aca_box_location_expressed = rules.h_aca_location.output.h_aca_box_location_expressed,
        h_aca_box_location_not_expressed = rules.h_aca_location.output.h_aca_box_location_not_expressed
    output:
        c_expressed = config['path']['c_box_expressed_fa'],
        c_not_expressed = config['path']['c_box_not_expressed_fa'],
        d_expressed = config['path']['d_box_expressed_fa'],
        d_not_expressed = config['path']['d_box_not_expressed_fa'],
        c_prime_expressed = config['path']['c_prime_box_expressed_fa'],
        c_prime_not_expressed = config['path']['c_prime_box_not_expressed_fa'],
        d_prime_expressed = config['path']['d_prime_box_expressed_fa'],
        d_prime_not_expressed = config['path']['d_prime_box_not_expressed_fa'],
        h_expressed = config['path']['h_box_expressed_fa'],
        h_not_expressed = config['path']['h_box_not_expressed_fa'],
        aca_expressed = config['path']['aca_box_expressed_fa'],
        aca_not_expressed = config['path']['aca_box_not_expressed_fa']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/convert_motif_table_to_fa.py"

rule logomaker:
    """ Create logo (consensus sequence) of c,d,c',d' boxes from sequences
        of expressed or not expressed C/D box snoRNAs. Same for H and ACA boxes."""
    input:
        box_fastas = rules.convert_motif_table_to_fa.output
    output:
        box_logos = expand(os.path.join(config['figures']['logo'], '{box}.svg'),
                            box=['expressed_c_box', 'not_expressed_c_box',
                                'expressed_d_box', 'not_expressed_d_box',
                                'expressed_c_prime_box', 'not_expressed_c_prime_box',
                                'expressed_d_prime_box', 'not_expressed_d_prime_box',
                                'expressed_h_box', 'not_expressed_h_box',
                                'expressed_aca_box', 'not_expressed_aca_box'])
    conda:
        "../envs/logomaker.yaml"
    script:
        "../scripts/python/graphs/logo_box.py"

rule logomaker_w_flanking_nt:
    """ Create logo (consensus sequence) of c,d,c',d' boxes from sequences
        of expressed or not expressed C/D box snoRNAs and with 3 flanking nt
        each side. Same for H and ACA boxes."""
    input:
        box_fastas = rules.flanking_nt_to_motifs.output
    output:
        box_logos = expand(os.path.join(config['figures']['logo_w_flanking'], '{box}.svg'),
                            box=['expressed_c_box_w_flanking_nt', 'not_expressed_c_box_w_flanking_nt',
                                'expressed_d_box_w_flanking_nt', 'not_expressed_d_box_w_flanking_nt',
                                'expressed_c_prime_box_w_flanking_nt', 'not_expressed_c_prime_box_w_flanking_nt',
                                'expressed_d_prime_box_w_flanking_nt', 'not_expressed_d_prime_box_w_flanking_nt',
                                'expressed_h_box_w_flanking_nt', 'not_expressed_h_box_w_flanking_nt',
                                'expressed_aca_box_w_flanking_nt', 'not_expressed_aca_box_w_flanking_nt'])
    conda:
        "../envs/logomaker.yaml"
    script:
        "../scripts/python/graphs/logo_box.py"
