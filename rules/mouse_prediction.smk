import os

include: "tgirt_seq_mouse.smk"
include: "cv_train_test_manual_split_top3.smk"
include: "cv_train_test_manual_split_top4.smk"
include: "downloads.smk"

rule mouse_gtf_to_bed:
    """ Convert a gtf file into a bed file (first command) and generate out of
        it a bed file of all snoRNA genes (second command)."""
    input:
        gtf = config['path']['mouse_gtf']
    output:
        gtf_bed = config['path']['sorted_gtf_bed_mouse'],
        all_sno_bed = config['path']['all_sno_bed_mouse']
    shell:
        """awk -v OFS='\t' 'NR>6 {{print $1, $4, $5, "to_remove"$10"to_remove", $6, $7, $2, $3, $8, "to_delete"$0}}' {input.gtf} | """
        """sed -E 's/to_remove"//g; s/";to_remove//g; s/to_delete.*gene_id/gene_id/g' | """
        """sort -n -k1,1 -k2,2 > {output.gtf_bed} && """
        """awk '$8=="gene" {{print $0}}' {output.gtf_bed}  | grep snoRNA | sed 's/\t$//g; s/^/chr/g' | sort -k1,1 -k2,2n > {output.all_sno_bed}"""


rule fasta_sno_sequence_mouse:
    """ Create a fasta of all mouse snoRNA sequences."""
    input:
        genome = rules.ensembl_mouse_genome.output.genome,
        sno_bed = rules.mouse_gtf_to_bed.output.all_sno_bed
    output:
        sno_fasta = config['path']['sno_sequences_mouse']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/fasta_sno_sequence_species.py"


rule find_mouse_snoRNA_type:
    """ Based on RNAcentral information, determine snoRNA type of all snoRNAs in
        mouse. Return this table merged to the tpm df computed by TGIRT-Seq."""
    input:
        rna_central_df = rules.get_RNA_central_snoRNAs.output.RNA_central_snoRNAs,
        id_conversion_df = rules.rna_central_to_ensembl_id.output.conversion_table,
        tpm_df = rules.merge_coco_output_mouse.output.merged_tpm,
        gtf = config['path']['mouse_gtf'],
        sno_fasta = rules.fasta_sno_sequence_mouse.output.sno_fasta
    output:
        snoRNA_type_df = config['path']['mouse_snoRNA_type']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/find_mouse_snoRNA_type.py"


rule find_mouse_snoRNA_labels_w_length:
    """ From the TGIRT-Seq pipeline TPM output, define snoRNAs as expressed or
        not expressed. Gather also the snoRNA length."""
    input:
        sno_tpm_df = rules.find_mouse_snoRNA_type.output.snoRNA_type_df,
        sno_bed = rules.mouse_gtf_to_bed.output.all_sno_bed
    output:
        tpm_label_df = config['path']['mouse_snoRNA_labels_w_length']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/find_mouse_snoRNA_labels_w_length.py"

rule format_gtf_bed_for_HG:
    """ Keep only gene features and remove snoRNA and other embedded genes
        from a gtf converted to bed format."""
    input:
        gtf_bed = rules.mouse_gtf_to_bed.output.gtf_bed
    output:
        formatted_gtf_bed = config['path']['formatted_gtf_bed_mouse']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/format_gtf_bed_for_HG.py"

rule find_mouse_snoRNA_HG:
    """ Intersect the snoRNA bed file and the gtf to see which snoRNAs are
        intronic and which are intergenic. Then return the host gene id."""
    input:
        formatted_gtf_bed = rules.format_gtf_bed_for_HG.output.formatted_gtf_bed,
        sno_bed = rules.mouse_gtf_to_bed.output.all_sno_bed
    output:
        mouse_snoRNA_HG = config['path']['mouse_snoRNA_HG']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/find_mouse_snoRNA_HG.py"


rule find_mouse_HG_expression_level:
    """ Get the abundance_cutoff for mouse snoRNA host genes. The RNA-Seq
        samples are from Shen et al. 2012 Nature and processed using the
        Recount3 pipeline (Wilks et al. 2021 Genome Biology)."""
    input:
        host_df = rules.find_mouse_snoRNA_HG.output.mouse_snoRNA_HG,
        sno_tpm_df = rules.find_mouse_snoRNA_labels_w_length.output.tpm_label_df,
        recount_tpm_df = rules.dowload_mouse_HG_RNA_seq_datasets.output.dataset
    output:
        HG_abundance_df = config['path']['mouse_HG_abundance_cutoff']
    params:
        srr_id_conversion = config['srr_id_conversion']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/find_mouse_HG_expression_level.py"


rule structure_mouse:
    """ Get the structure stability of mouse snoRNAs (their MFE or Minimal Free
        Energy) using RNAfold (first convert T to U in the sequences). """
    input:
        sequences = rules.fasta_sno_sequence_mouse.output.sno_fasta
    params:
        temp_name = "sno_mouse.mfe"
    output:
        mfe = config['path']['structure_stability_fasta_mouse']
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "sed -E '/^[^>]/ s/T/U/g' {input.sequences} | sed 's/>/>MOUSE_/g' > temp_structure && "
        "RNAfold --infile=temp_structure --outfile={params.temp_name} && sed -i 's/MOUSE_//g' {params.temp_name} && "
        "mv {params.temp_name} {output.mfe} && mkdir -p data/structure/stability_mouse/ && "
        "mv MOUSE*.ps data/structure/stability_mouse/ && rm temp_structure"

rule fasta_to_tsv_mouse:
    """ Convert the fasta output of RNA fold into a tsv table with a snoRNA id
        as the first column and the MFE as the second column."""
    input:
        mfe = rules.structure_mouse.output.mfe
    params:
        temp_mfe = "temporary_mouse_mfe.tsv",
        temp_id = "temporary_mouse_ids.tsv"
    output:
        mfe_final = config['path']['structure_stability_tsv_mouse']
    shell:
        """grep -E ">" {input.mfe} | sed 's/>//g' > {params.temp_id} && """
        """grep -oE "\-*[0-9]+\.[0-9]*" {input.mfe} > {params.temp_mfe} && """
        """paste {params.temp_id} {params.temp_mfe} > {output.mfe_final} && """
        """rm temporary_mouse*"""

rule get_chr_size:
    """ Get the .genome file containing the size of each chromosome using
        samtools faidx and the fasta of the genome. Samtools faidx creates a
        .fai file of the genome containing the chr size in the 2nd column."""
    input:
        genome = rules.ensembl_mouse_genome.output.genome
    output:
        genome_chr_size = 'data/references/mouse_genome_chr_sizes.genome'
    conda:
        "../envs/samtools.yaml"
    shell:
        "mkdir -p log/ && samtools faidx {input.genome} && "
        "cut -f1,2 {input.genome}.fai > {output.genome_chr_size}"

rule flank_extend_snoRNA_mouse:
    """ Get the upstream and downstream flanking regions (15 nt) of each snoRNA
        using pybedtools flank. Then extend these flanking regions inside the
        snoRNA sequence using pybedtools slop. Extend for 5 nt from the 5' and
        3' of C/D box snoRNAs; extend 5 nt from the 5' and 3 nt from the 3' of
        H/ACA box snoRNAs."""
    input:
        all_sno_bed = rules.mouse_gtf_to_bed.output.all_sno_bed,
        sno_info = rules.find_mouse_snoRNA_type.output.snoRNA_type_df,
        genome_chr_size = rules.get_chr_size.output.genome_chr_size
    output:
        flanking_cd_left = config['path']['flanking_regions_cd_left_mouse'],
        flanking_cd_right = config['path']['flanking_regions_cd_right_mouse'],
        flanking_haca_left = config['path']['flanking_regions_haca_left_mouse'],
        flanking_haca_right = config['path']['flanking_regions_haca_right_mouse']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/flank_extend_snoRNA_mouse.py"

rule get_fasta_terminal_stem_mouse:
    """ Get the sequence of all the extended flanking regions into a fasta file
        using pybedtools 'sequence'."""
    input:
        genome_fasta = rules.ensembl_mouse_genome.output.genome,
        flanking_cd_left = rules.flank_extend_snoRNA_mouse.output[0],
        flanking_cd_right = rules.flank_extend_snoRNA_mouse.output[1],
        flanking_haca_left = rules.flank_extend_snoRNA_mouse.output[2],
        flanking_haca_right = rules.flank_extend_snoRNA_mouse.output[3]
    output:
        sequences = config['path']['flanking_regions_sno_fasta_mouse']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_fasta_terminal_stem.py"


rule rna_cofold_mouse:
    """ Get the structure stability of the potential terminal stem of snoRNAs."""
    input:
        fasta = rules.get_fasta_terminal_stem_mouse.output.sequences
    output:
        mfe_stem = config['path']['terminal_stem_mfe_fasta_mouse']
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "sed 's/>/>Mouse_/g' {input.fasta} > Mouse_cofold.fa && "
        "RNAcofold < Mouse_cofold.fa > {output.mfe_stem} && sed -i 's/Mouse_//g' {output.mfe_stem} && "
        "mkdir -p data/terminal_stem_mouse/ && mv Mouse*.ps data/terminal_stem_mouse/ && rm Mouse_cofold.fa"


rule fasta_to_tsv_terminal_stem_mfe_mouse:
    """ Convert the fasta output of RNAcofold into a tsv table with a snoRNA id
        as the first column and the MFE of the terminal stem as the second column."""
    input:
        mfe = rules.rna_cofold_mouse.output.mfe_stem
    params:
        temp_mfe = "temporary_terminal_stem_mfe.tsv",
        temp_id = "temporary_terminal_stem_mfe_ids.tsv"
    output:
        mfe_stem_final = config['path']['terminal_stem_mfe_tsv_mouse']
    shell:
        """grep -E ">" {input.mfe} | sed 's/>//g' > {params.temp_id} && """
        """grep -oE "\-*[0-9]+\.[0-9]*" {input.mfe} > {params.temp_mfe} && """
        """paste {params.temp_id} {params.temp_mfe} > {output.mfe_stem_final} && """
        """rm temporary_terminal_stem*"""


rule fasta_per_sno_type_mouse:
    """ Create a fasta per snoRNA type out of of a fasta containing all snoRNAs."""
    input:
        sno_fasta = rules.fasta_sno_sequence_mouse.output.sno_fasta,
        sno_info = rules.find_mouse_snoRNA_type.output.snoRNA_type_df
    output:
        cd_fasta = config['path']['all_cd_mouse_fa'],
        haca_fasta = config['path']['all_haca_mouse_fa']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/fasta_per_sno_type_mouse.py"


rule c_d_box_location_all_mouse:
    """ Find C, D, C' and D' motif and location (if they exist) in ALL C/D box
        snoRNAs. The output of this rule will be used to create the hamming
        distance feature for all boxes of C/D snoRNAs in the initial feature_df."""
    input:
        cd_fasta = rules.fasta_per_sno_type_mouse.output.cd_fasta
    output:
        c_d_box_location = config['path']['c_d_and_prime_box_location_mouse']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/cd_box_location_all.py"


rule h_aca_box_location_all_mouse:
    """ Find H and ACA motif and location (if they exist) in ALL H/ACA box
        snoRNAs. The output of this rule will be used to create the hamming
        distance feature for all boxes of H/ACA snoRNAs in the initial feature_df."""
    input:
        dot_bracket = rules.structure_mouse.output.mfe,
        haca_fasta = rules.fasta_per_sno_type_mouse.output.haca_fasta
    output:
        h_aca_box_location = config['path']['h_aca_box_location_mouse']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/haca_box_location_all.py"

rule hamming_distance_box_all_mouse:
    """ Get the hamming distance per box and also a combined hamming distance
        per snoRNA boxes for all snoRNAs as input features of the models."""
    input:
        c_d_box_location = rules.c_d_box_location_all_mouse.output.c_d_box_location,
        h_aca_box_location = rules.h_aca_box_location_all_mouse.output.h_aca_box_location
    output:
        hamming_distance_box_df = config['path']['hamming_distance_box_df_mouse']
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/hamming_distance_box_all.py"


rule merge_features_label_mouse:
    """ Merge all feature columns and label inside one dataframe."""
    input:
        abundance_cutoff = rules.find_mouse_HG_expression_level.output.HG_abundance_df,
        sno_structure_mfe = rules.fasta_to_tsv_mouse.output.mfe_final,
        terminal_stem_mfe = rules.fasta_to_tsv_terminal_stem_mfe_mouse.output.mfe_stem_final,
        hamming_distance_box = rules.hamming_distance_box_all_mouse.output.hamming_distance_box_df
    output:
        feature_df = config['path']['feature_df_mouse']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/merge_features_mouse.py"


rule predict_mouse_snoRNA_label:
    """ Predict the abundance status of all mouse snoRNA based on the
        top4 features (combined_box_hamming. sno_mfe, terminal_stem_mfe and
        host_expressed). We scale the feature_df using the same parameters (mean
        and stdev) used to scale each training set."""
    input:
        feature_df = rules.merge_features_label_mouse.output.feature_df,
        human_snoRNA_feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        pickled_trained_model = rules.train_models_scale_after_manual_split_top4.output.pickled_trained_model
    output:
        predicted_label_df = 'results/tables/mouse_prediction/{models2}_predicted_label_{manual_iteration}.tsv',
        scaled_feature_df = 'results/tables/mouse_prediction/{models2}_scaled_features_{manual_iteration}.tsv',
        label_df = 'results/tables/mouse_prediction/{models2}_label_{manual_iteration}.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/predict_mouse_snoRNA_label.py"

rule test_accuracy_mouse:
    """ Test model performance on mouse snoRNA data and return their accuracies."""
    input:
        X_test = rules.predict_mouse_snoRNA_label.output.scaled_feature_df,
        y_test = rules.predict_mouse_snoRNA_label.output.label_df,
        pickled_trained_model = rules.train_models_scale_after_manual_split_top4.output.pickled_trained_model
    output:
        test_accuracy = os.path.join(config['path']['test_accuracy_mouse'],
                                    '{models2}_test_accuracy_{manual_iteration}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/test_models_scale_after_split.py"

rule confusion_matrix_f1_mouse:
    """ Create a confusion matrix (all True/False positives or negatives) for
        each model in order to compare whether they misclassify the same snoRNAs
        or not. Compute also the F1 score per model. Do this on mouse snoRNAs."""
    input:
        X_test = rules.predict_mouse_snoRNA_label.output.scaled_feature_df,
        y_test = rules.predict_mouse_snoRNA_label.output.label_df,
        pickled_trained_model = rules.train_models_scale_after_manual_split_top4.output.pickled_trained_model
    output:
        confusion_matrix = os.path.join(config['path']['confusion_matrix_f1_mouse'],
                        '{models2}_confusion_matrix_w_f1_score_{manual_iteration}.tsv'),
        info_df = os.path.join(config['path']['confusion_matrix_f1_mouse'],
                        '{models2}_confusion_matrix_values_per_sno_{manual_iteration}.tsv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/confusion_matrix_f1_scale_after_split.py"
