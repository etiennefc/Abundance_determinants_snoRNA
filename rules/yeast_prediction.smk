import os

include: "cv_train_test_manual_split_top3.smk"
include: "cv_train_test_manual_split_top4.smk"
include: "downloads.smk"

rule yeast_gtf_to_bed:
    """ Convert a gtf file into a bed file (first command) and generate out of
        it a bed file of all snoRNA genes (second command)."""
    input:
        gtf = rules.download_yeast_annotation.output.std_gtf
    output:
        gtf_bed = config['path']['sorted_gtf_bed_yeast'],
        all_sno_bed = config['path']['all_sno_bed_yeast']
    shell:
        """awk -v OFS='\t' 'NR>6 {{print $1, $4, $5, "to_remove"$10"to_remove", $6, $7, $2, $3, $8, "to_delete"$0}}' {input.gtf} | """
        """sed -E 's/to_remove"//g; s/";to_remove//g; s/to_delete.*gene_id/gene_id/g' | """
        """sort -n -k1,1 -k2,2 > {output.gtf_bed} && """
        """awk '$8=="gene" {{print $0}}' {output.gtf_bed}  | grep snoRNA | sed 's/\t$//g; s/^/chr/g' | sort -k1,1 -k2,2n > {output.all_sno_bed}"""


rule fasta_sno_sequence_yeast:
    """ Create a fasta of all yeast snoRNA sequences."""
    input:
        genome = rules.download_yeast_genome.output.genome,
        sno_bed = rules.yeast_gtf_to_bed.output.all_sno_bed
    output:
        sno_fasta = config['path']['sno_sequences_yeast']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/fasta_sno_sequence_species.py"


rule find_yeast_snoRNA_type:
    """ Based on Yeast mine information, determine snoRNA type of all snoRNAs in
        yeast. Return this table merged to the tpm df computed by TGIRT-Seq in 3
        WT yeast samples."""
    input:
        yeast_mine_df = config['path']['yeast_snoRNA_type'],
        tpm_df = 'data/references/merged_tpm_WT_yeast.tsv'
    output:
        snoRNA_type_df = 'results/tables/yeast_snoRNA_type_df.tsv'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/find_yeast_snoRNA_type.py"


rule find_yeast_snoRNA_labels_w_length:
    """ From the TGIRT-Seq pipeline TPM output, define snoRNAs as expressed or
        not expressed. Gather also the snoRNA length."""
    input:
        sno_tpm_df = rules.find_yeast_snoRNA_type.output.snoRNA_type_df,
        sno_bed = rules.yeast_gtf_to_bed.output.all_sno_bed
    output:
        tpm_label_df = config['path']['yeast_snoRNA_labels_w_length']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/find_mouse_snoRNA_labels_w_length.py"

rule format_gtf_bed_for_HG_yeast:
    """ Keep only gene features and remove snoRNA and other embedded genes
        from a gtf converted to bed format."""
    input:
        gtf_bed = rules.yeast_gtf_to_bed.output.gtf_bed
    output:
        formatted_gtf_bed = config['path']['formatted_gtf_bed_yeast']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/format_gtf_bed_for_HG.py"

rule find_yeast_snoRNA_HG:
    """ Intersect the snoRNA bed file and the gtf to see which snoRNAs are
        intronic and which are intergenic. Then return the host gene id."""
    input:
        formatted_gtf_bed = rules.format_gtf_bed_for_HG_yeast.output.formatted_gtf_bed,
        sno_bed = rules.yeast_gtf_to_bed.output.all_sno_bed
    output:
        species_snoRNA_HG = config['path']['yeast_snoRNA_HG']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/find_species_snoRNA_HG.py"


rule find_yeast_HG_expression_level:
    """ Get the abundance_cutoff for yeast snoRNA host genes. The TGIRT-Seq data
        is from 3 wild type S. cerevisiae samples."""
    input:
        host_df = rules.find_yeast_snoRNA_HG.output.species_snoRNA_HG,
        sno_tpm_df = rules.find_yeast_snoRNA_labels_w_length.output.tpm_label_df,
        tpm_df = 'data/references/merged_tpm_WT_yeast.tsv'
    output:
        HG_abundance_df = config['path']['yeast_HG_abundance_cutoff']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/find_yeast_HG_expression_level.py"

rule structure_yeast:
    """ Get the structure stability of yeast snoRNAs (their MFE or Minimal Free
        Energy) using RNAfold (first convert T to U in the sequences). """
    input:
        sequences = rules.fasta_sno_sequence_yeast.output.sno_fasta
    params:
        temp_name = "sno_yeast.mfe"
    output:
        mfe = config['path']['structure_stability_fasta_yeast']
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "mkdir -p data/structure/stability_yeast/ && "
        "sed -E '/^[^>]/ s/T/U/g' {input.sequences} > temp_structure_yeast && "
        "RNAfold --infile=temp_structure_yeast --outfile={params.temp_name} && "
        "mv {params.temp_name} {output.mfe} && "
        "mv *.ps data/structure/stability_yeast/ && rm temp_structure_yeast"

rule fasta_to_tsv_yeast:
    """ Convert the fasta output of RNA fold into a tsv table with a snoRNA id
        as the first column and the MFE as the second column."""
    input:
        mfe = rules.structure_yeast.output.mfe
    params:
        temp_mfe = "temporary_yeast_mfe.tsv",
        temp_id = "temporary_yeast_ids.tsv"
    output:
        mfe_final = config['path']['structure_stability_tsv_yeast']
    shell:
        """grep -E ">" {input.mfe} | sed 's/>//g' > {params.temp_id} && """
        """grep -oE "\-*[0-9]+\.[0-9]*" {input.mfe} > {params.temp_mfe} && """
        """paste {params.temp_id} {params.temp_mfe} > {output.mfe_final} && """
        """rm temporary_yeast*"""

rule get_chr_size_yeast:
    """ Get the .genome file containing the size of each chromosome using
        samtools faidx and the fasta of the genome. Samtools faidx creates a
        .fai file of the genome containing the chr size in the 2nd column."""
    input:
        genome = rules.download_yeast_genome.output.genome
    output:
        genome_chr_size = 'data/references/yeast_genome_chr_sizes.genome'
    conda:
        "../envs/samtools.yaml"
    shell:
        "mkdir -p log/ && samtools faidx {input.genome} && "
        "cut -f1,2 {input.genome}.fai > {output.genome_chr_size}"

rule flank_extend_snoRNA_yeast:
    """ Get the upstream and downstream flanking regions (15 nt) of each snoRNA
        using pybedtools flank. Then extend these flanking regions inside the
        snoRNA sequence using pybedtools slop. Extend for 5 nt from the 5' and
        3' of C/D box snoRNAs; extend 5 nt from the 5' and 3 nt from the 3' of
        H/ACA box snoRNAs."""
    input:
        all_sno_bed = rules.yeast_gtf_to_bed.output.all_sno_bed,
        sno_info = rules.find_yeast_snoRNA_type.output.snoRNA_type_df,
        genome_chr_size = rules.get_chr_size_yeast.output.genome_chr_size
    output:
        flanking_cd_left = config['path']['flanking_regions_cd_left_yeast'],
        flanking_cd_right = config['path']['flanking_regions_cd_right_yeast'],
        flanking_haca_left = config['path']['flanking_regions_haca_left_yeast'],
        flanking_haca_right = config['path']['flanking_regions_haca_right_yeast']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/flank_extend_snoRNA_mouse.py"

rule get_fasta_terminal_stem_yeast:
    """ Get the sequence of all the extended flanking regions into a fasta file
        using pybedtools 'sequence'."""
    input:
        genome_fasta = rules.download_yeast_genome.output.genome,
        flanking_cd_left = rules.flank_extend_snoRNA_yeast.output[0],
        flanking_cd_right = rules.flank_extend_snoRNA_yeast.output[1],
        flanking_haca_left = rules.flank_extend_snoRNA_yeast.output[2],
        flanking_haca_right = rules.flank_extend_snoRNA_yeast.output[3]
    output:
        sequences = config['path']['flanking_regions_sno_fasta_yeast']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_fasta_terminal_stem.py"


rule rna_cofold_yeast:
    """ Get the structure stability of the potential terminal stem of snoRNAs."""
    input:
        fasta = rules.get_fasta_terminal_stem_yeast.output.sequences
    output:
        mfe_stem = config['path']['terminal_stem_mfe_fasta_yeast']
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "RNAcofold < {input.fasta} > {output.mfe_stem} && "
        "mv *.ps data/terminal_stem_yeast/"


rule fasta_to_tsv_terminal_stem_mfe_yeast:
    """ Convert the fasta output of RNAcofold into a tsv table with a snoRNA id
        as the first column and the MFE of the terminal stem as the second column."""
    input:
        mfe = rules.rna_cofold_yeast.output.mfe_stem
    params:
        temp_mfe = "temporary_terminal_stem_mfe_yeast.tsv",
        temp_id = "temporary_terminal_stem_mfe_yeast_ids.tsv"
    output:
        mfe_stem_final = config['path']['terminal_stem_mfe_tsv_yeast']
    shell:
        """grep -E ">" {input.mfe} | sed 's/>//g' > {params.temp_id} && """
        """grep -oE "\-*[0-9]+\.[0-9]*" {input.mfe} > {params.temp_mfe} && """
        """paste {params.temp_id} {params.temp_mfe} > {output.mfe_stem_final} && """
        """rm temporary_terminal_stem_mfe_yeast*"""


rule fasta_per_sno_type_yeast:
    """ Create a fasta per snoRNA type out of of a fasta containing all snoRNAs."""
    input:
        sno_fasta = rules.fasta_sno_sequence_yeast.output.sno_fasta,
        sno_info = rules.find_yeast_snoRNA_type.output.snoRNA_type_df
    output:
        cd_fasta = config['path']['all_cd_yeast_fa'],
        haca_fasta = config['path']['all_haca_yeast_fa']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/fasta_per_sno_type_mouse.py"


rule c_d_box_location_all_yeast:
    """ Find C, D, C' and D' motif and location (if they exist) in ALL C/D box
        snoRNAs. The output of this rule will be used to create the hamming
        distance feature for all boxes of C/D snoRNAs in the initial feature_df."""
    input:
        cd_fasta = rules.fasta_per_sno_type_yeast.output.cd_fasta
    output:
        c_d_box_location = config['path']['c_d_and_prime_box_location_yeast']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/cd_box_location_all.py"


rule h_aca_box_location_all_yeast:
    """ Find H and ACA motif and location (if they exist) in ALL H/ACA box
        snoRNAs. The output of this rule will be used to create the hamming
        distance feature for all boxes of H/ACA snoRNAs in the initial feature_df."""
    input:
        dot_bracket = rules.structure_yeast.output.mfe,
        haca_fasta = rules.fasta_per_sno_type_yeast.output.haca_fasta
    output:
        h_aca_box_location = config['path']['h_aca_box_location_yeast']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/haca_box_location_all.py"

rule hamming_distance_box_all_yeast:
    """ Get the hamming distance per box and also a combined hamming distance
        per snoRNA boxes for all snoRNAs as input features of the models."""
    input:
        c_d_box_location = rules.c_d_box_location_all_yeast.output.c_d_box_location,
        h_aca_box_location = rules.h_aca_box_location_all_yeast.output.h_aca_box_location
    output:
        hamming_distance_box_df = config['path']['hamming_distance_box_df_yeast']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/hamming_distance_box_all.py"


rule merge_features_label_yeast:
    """ Merge all feature columns and label inside one dataframe."""
    input:
        abundance_cutoff = rules.find_yeast_HG_expression_level.output.HG_abundance_df,
        sno_structure_mfe = rules.fasta_to_tsv_yeast.output.mfe_final,
        terminal_stem_mfe = rules.fasta_to_tsv_terminal_stem_mfe_yeast.output.mfe_stem_final,
        hamming_distance_box = rules.hamming_distance_box_all_yeast.output.hamming_distance_box_df
    output:
        feature_df = config['path']['feature_df_yeast']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/merge_features_mouse.py"


rule predict_yeast_snoRNA_label:
    """ Predict the abundance status of yeast snoRNA based on the
        top4 features (combined_box_hamming, sno_mfe, terminal_stem_mfe and
        host_expressed). We scale the feature_df using the same parameters (mean
        and stdev) used to scale each training set in human. We use the same log_reg
        threshold that was defined on mouse predictions to maximize true
        positive rate and minimize the false positive rate. Return the predicted
        label and all features columns."""
    input:
        feature_df = rules.merge_features_label_yeast.output.feature_df,
        human_snoRNA_feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        pickled_trained_model = expand(rules.train_test_accuracy_species_prediction_top4_log_reg_thresh.output.pickled_trained_model, rs=42),
        threshold = expand(rules.train_test_accuracy_species_prediction_top4_log_reg_thresh.output.threshold, rs=42)
    output:
        predicted_label_df = 'results/tables/yeast_prediction/yeast_predicted_label.tsv',
        scaled_feature_df = 'results/tables/yeast_prediction/yeast_scaled_features.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/predict_species_snoRNA_label_final.py"

rule predict_yeast_snoRNA_label_no_thresh:
    """ Predict the abundance status of yeast snoRNA based on the
        top4 features (combined_box_hamming, sno_mfe, terminal_stem_mfe and
        host_expressed). We scale the feature_df using the same parameters (mean
        and stdev) used to scale each training set in human. We use the
        LogisticRegression model (hyperparameter_tuning and training on 10% and
        90% of human snoRNAs respectively). Return the predicted
        label and all features columns."""
    input:
        feature_df = rules.merge_features_label_yeast.output.feature_df,
        human_snoRNA_feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        pickled_trained_model = expand(rules.train_models_species_prediction_top4_rs.output.pickled_trained_model, rs=42, models2='log_reg')
    output:
        predicted_label_df = 'results/tables/yeast_prediction/yeast_predicted_label_no_thresh.tsv',
        scaled_feature_df = 'results/tables/yeast_prediction/yeast_scaled_features_no_thresh.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/predict_yeast_snoRNA_label.py"
