import os

include: "downloads_species.smk"

rule species_gtf_to_bed:
    """ Convert a gtf file into a bed file (first command) and generate out of
        it a bed file of all snoRNA genes (second command)."""
    input:
        gtf = rules.ensembl_species_gtf.output.gtf
    output:
        gtf_bed = 'data/bed_files/{species}_sorted.gtf.bed',
        all_sno_bed = 'data/bed_files/all_{species}_snoRNA.bed'
    shell:
        """awk -v OFS='\t' 'NR>6 {{print $1, $4, $5, "to_remove"$10"to_remove", $6, $7, $2, $3, $8, "to_delete"$0}}' {input.gtf} | """
        """sed -E 's/to_remove"//g; s/";to_remove//g; s/to_delete.*gene_id/gene_id/g' | """
        """sort -n -k1,1 -k2,2 > {output.gtf_bed} && """
        """awk '$8=="gene" {{print $0}}' {output.gtf_bed}  | grep snoRNA | sed 's/\t$//g; s/^/chr/g' | sort -k1,1 -k2,2n > {output.all_sno_bed}"""

rule get_chr_size_species:
    """ Get the .genome file containing the size of each chromosome using
        samtools faidx and the fasta of the genome. Samtools faidx creates a
        .fai file of the genome containing the chr size in the 2nd column."""
    input:
        genome = rules.ensembl_species_genome.output.genome
    output:
        genome_chr_size = 'data/references/{species}_genome_chr_sizes.genome'
    conda:
        "../envs/samtools.yaml"
    shell:
        "mkdir -p log/ && samtools faidx {input.genome} && "
        "cut -f1,2 {input.genome}.fai > {output.genome_chr_size}"

rule filter_genome:
    """ Filter out weird chromosomes from fasta of genome (especially for
        haploid (ALT) chromosomes in Danio rerio)."""
    input:
        genome = rules.ensembl_species_genome.output.genome,
        chr_size = rules.get_chr_size_species.output.genome_chr_size
    output:
        filtered_genome = 'data/references/{species}_filtered_genome.fa'
    params:
        temp_filter = 'temp_filter_{species}'
    conda:
        "../envs/samtools.yaml"
    shell:
        """for i in $(cut -f1 {input.chr_size}); do if [[ $i != *"_ALT_"* ]]; """
        """then echo $i && samtools faidx {input.genome} $i >> {params.temp_filter}; fi; done && """
        """sed 's/>/>chr/' {params.temp_filter} > {output.filtered_genome} && rm {params.temp_filter}"""

rule get_chr_size_filtered_species:
    """ Get the .genome file containing the size of each filtered chromosome using
        samtools faidx and the fasta of the genome. Samtools faidx creates a
        .fai file of the genome containing the chr size in the 2nd column."""
    input:
        genome = rules.filter_genome.output.filtered_genome
    output:
        genome_chr_size = 'data/references/{species}_genome_filtered_chr_sizes.genome'
    conda:
        "../envs/samtools.yaml"
    shell:
        "mkdir -p log/ && samtools faidx {input.genome} && "
        "cut -f1,2 {input.genome}.fai > {output.genome_chr_size}"

rule fasta_sno_sequence_species:
    """ Create a fasta of all species snoRNA sequences."""
    input:
        genome = rules.filter_genome.output.filtered_genome,
        sno_bed = rules.species_gtf_to_bed.output.all_sno_bed,
        fake_dependency = rules.get_chr_size_filtered_species.output.genome_chr_size
    output:
        sno_fasta = 'data/bed_files/sno_sequences_{species}.fa'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/fasta_sno_sequence_species.py"

rule find_species_snoRNA_type:
    """ Based on RNAcentral information, determine snoRNA type of all snoRNAs in
        species Ensembl gtf."""
    input:
        rna_central_df = rules.get_RNA_central_snoRNAs_species.output.RNA_central_snoRNAs,
        id_conversion_df = rules.rna_central_to_ensembl_id_species.output.conversion_table,
        sno_bed = rules.species_gtf_to_bed.output.all_sno_bed
    output:
        snoRNA_type_df = "results/tables/{species}_snoRNA_type.tsv"
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/find_species_snoRNA_type.py"

rule format_gtf_bed_for_HG_species:
    """ Keep only gene features and remove snoRNA and other embedded genes
        from a gtf converted to bed format."""
    input:
        gtf_bed = rules.species_gtf_to_bed.output.gtf_bed
    output:
        formatted_gtf_bed = "data/bed_files/{species}_formatted.gtf.bed"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/format_gtf_bed_for_HG.py"

rule find_species_snoRNA_HG:
    """ Intersect the snoRNA bed file and the gtf to see which snoRNAs are
        intronic and which are intergenic. Then return the host gene id."""
    input:
        formatted_gtf_bed = rules.format_gtf_bed_for_HG_species.output.formatted_gtf_bed,
        sno_bed = rules.species_gtf_to_bed.output.all_sno_bed
    output:
        species_snoRNA_HG = "data/references/{species}_host_gene_list.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/find_species_snoRNA_HG.py"

rule find_species_HG_expression_level:
    """ Get the abundance_cutoff for species snoRNA host genes. The RNA-Seq
        samples are from the Bgee database of normal animal samples (Bastian et
        al., NAR, 2021)."""
    input:
        host_df = rules.find_species_snoRNA_HG.output.species_snoRNA_HG,
        sno_df = rules.find_species_snoRNA_type.output.snoRNA_type_df,
        tpm_df_dir = rules.get_species_HG_expression.output.expression_dir
    output:
        HG_abundance_df = "data/references/{species}_HG_abundance_cutoff.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/find_species_HG_expression_level.py"

rule structure_species:
    """ Get the structure stability of species snoRNAs (their MFE or Minimal Free
        Energy) using RNAfold (first convert T to U in the sequences). """
    input:
        sequences = rules.fasta_sno_sequence_species.output.sno_fasta
    params:
        temp_name = "sno_{species}.mfe",
        mfe_dir = "data/structure/stability_{species}"
    output:
        mfe = "data/structure/stability_{species}/sno_mfe_{species}.fa"
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "sed -E '/^[^>]/ s/T/U/g' {input.sequences} > {params.mfe_dir}/temp_structure_species && "
        "RNAfold --infile={params.mfe_dir}/temp_structure_species --outfile={params.temp_name} && "
        "mkdir -p {params.mfe_dir} && mv {params.temp_name} {output.mfe} && "
        "mv *.ps {params.mfe_dir} && rm {params.mfe_dir}/temp_structure_species"

rule fasta_to_tsv_species:
    """ Convert the fasta output of RNA fold into a tsv table with a snoRNA id
        as the first column and the MFE as the second column."""
    input:
        mfe = rules.structure_species.output.mfe
    params:
        temp_mfe = "temporary_{species}_mfe.tsv",
        temp_id = "temporary_{species}_ids.tsv"
    output:
        mfe_final = "data/structure/stability_{species}/sno_mfe_{species}.tsv"
    shell:
        """grep -E ">" {input.mfe} | sed 's/>//g' > {params.temp_id} && """
        """grep -oE "\-*[0-9]+\.[0-9]*" {input.mfe} > {params.temp_mfe} && """
        """paste {params.temp_id} {params.temp_mfe} > {output.mfe_final} && """
        """rm {params.temp_mfe} {params.temp_id}"""

rule flank_extend_snoRNA_species:
    """ Get the upstream and downstream flanking regions (15 nt) of each snoRNA
        using pybedtools flank. Then extend these flanking regions inside the
        snoRNA sequence using pybedtools slop. Extend for 5 nt from the 5' and
        3' of C/D box snoRNAs; extend 5 nt from the 5' and 3 nt from the 3' of
        H/ACA box snoRNAs."""
    input:
        all_sno_bed = rules.species_gtf_to_bed.output.all_sno_bed,
        sno_info = rules.find_species_snoRNA_type.output.snoRNA_type_df,
        genome_chr_size = rules.get_chr_size_filtered_species.output.genome_chr_size
    output:
        flanking_cd_left = "data/terminal_stem_{species}/extended_region_cd_left.bed",
        flanking_cd_right = "data/terminal_stem_{species}/extended_region_cd_right.bed",
        flanking_haca_left = "data/terminal_stem_{species}/extended_region_haca_left.bed",
        flanking_haca_right = "data/terminal_stem_{species}/extended_region_haca_right.bed"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/flank_extend_snoRNA_species.py"

rule get_fasta_terminal_stem_species:
    """ Get the sequence of all the extended flanking regions into a fasta file
        using pybedtools 'sequence'."""
    input:
        genome_fasta = rules.filter_genome.output.filtered_genome,
        flanking_cd_left = rules.flank_extend_snoRNA_species.output[0],
        flanking_cd_right = rules.flank_extend_snoRNA_species.output[1],
        flanking_haca_left = rules.flank_extend_snoRNA_species.output[2],
        flanking_haca_right = rules.flank_extend_snoRNA_species.output[3]
    output:
        sequences = "data/terminal_stem_{species}/terminal_stem_sequences.fa"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_fasta_terminal_stem.py"


rule rna_cofold_species:
    """ Get the structure stability of the potential terminal stem of snoRNAs."""
    input:
        fasta = rules.get_fasta_terminal_stem_species.output.sequences
    output:
        mfe_stem = "data/terminal_stem_{species}/terminal_stem_mfe.fa"
    params:
        mfe_dir = "data/terminal_stem_{species}/"
    conda:
        "../envs/rna_fold.yaml"
    shell:
        "RNAcofold < {input.fasta} > {output.mfe_stem} && "
        "mv *.ps {params.mfe_dir}"


rule fasta_to_tsv_terminal_stem_mfe_species:
    """ Convert the fasta output of RNAcofold into a tsv table with a snoRNA id
        as the first column and the MFE of the terminal stem as the second column."""
    input:
        mfe = rules.rna_cofold_species.output.mfe_stem
    params:
        temp_mfe = "temporary_terminal_stem_mfe_{species}.tsv",
        temp_id = "temporary_terminal_stem_mfe_ids_{species}.tsv"
    output:
        mfe_stem_final = "data/terminal_stem_{species}/terminal_stem_mfe.tsv"
    shell:
        """grep -E ">" {input.mfe} | sed 's/>//g' > {params.temp_id} && """
        """grep -oE "\-*[0-9]+\.[0-9]*" {input.mfe} > {params.temp_mfe} && """
        """paste {params.temp_id} {params.temp_mfe} > {output.mfe_stem_final} && """
        """rm {params.temp_mfe} {params.temp_id}"""


rule fasta_per_sno_type_species:
    """ Create a fasta per snoRNA type out of a fasta containing all snoRNAs."""
    input:
        sno_fasta = rules.fasta_sno_sequence_species.output.sno_fasta,
        sno_info = rules.find_species_snoRNA_type.output.snoRNA_type_df
    output:
        cd_fasta = "data/references/all_cd_{species}.fa",
        haca_fasta = "data/references/all_haca_{species}.fa"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/fasta_per_sno_type_species.py"


rule c_d_box_location_all_species:
    """ Find C, D, C' and D' motif and location (if they exist) in ALL C/D box
        snoRNAs. The output of this rule will be used to create the hamming
        distance feature for all boxes of C/D snoRNAs in the initial feature_df."""
    input:
        cd_fasta = rules.fasta_per_sno_type_species.output.cd_fasta
    output:
        c_d_box_location = "data/references/c_d_and_prime_box_location_{species}.tsv"
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/cd_box_location_all.py"


rule h_aca_box_location_all_species:
    """ Find H and ACA motif and location (if they exist) in ALL H/ACA box
        snoRNAs. The output of this rule will be used to create the hamming
        distance feature for all boxes of H/ACA snoRNAs in the initial feature_df."""
    input:
        dot_bracket = rules.structure_species.output.mfe,
        haca_fasta = rules.fasta_per_sno_type_species.output.haca_fasta
    output:
        h_aca_box_location = "data/references/h_aca_box_location_{species}.tsv"
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/haca_box_location_all.py"

rule hamming_distance_box_all_species:
    """ Get the hamming distance per box and also a combined hamming distance
        per snoRNA boxes for all snoRNAs as input features of the models."""
    input:
        c_d_box_location = rules.c_d_box_location_all_species.output.c_d_box_location,
        h_aca_box_location = rules.h_aca_box_location_all_species.output.h_aca_box_location
    output:
        hamming_distance_box_df = "data/references/hamming_distance_box_df_{species}.tsv"
    conda:
        "../envs/python_regex.yaml"
    script:
        "../scripts/python/hamming_distance_box_all.py"


rule merge_features_species:
    """ Merge all feature columns inside one dataframe."""
    input:
        abundance_cutoff = rules.find_species_HG_expression_level.output.HG_abundance_df,
        sno_structure_mfe = rules.fasta_to_tsv_species.output.mfe_final,
        terminal_stem_mfe = rules.fasta_to_tsv_terminal_stem_mfe_species.output.mfe_stem_final,
        hamming_distance_box = rules.hamming_distance_box_all_species.output.hamming_distance_box_df
    output:
        feature_df = "results/tables/all_features_{species}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/merge_features_species.py"

rule predict_species_snoRNA_label:
    """ Predict the abundance status of all species snoRNA based on the
        top4 features (combined_box_hamming, sno_mfe, terminal_stem_mfe and
        host_expressed). We scale the feature_df using the same parameters (mean
        and stdev) used to scale each training set. We use the same log_reg
        threshold that was defined on mouse predictions to maximize true
        positive rate and minimize the false positive rate. Return the predicted
        label and all features columns."""
    input:
        feature_df = rules.merge_features_species.output.feature_df,
        human_snoRNA_feature_df = rules.one_hot_encode_before_split.output.one_hot_encoded_df,
        pickled_trained_model = expand(rules.train_test_accuracy_species_prediction_top4_log_reg_thresh.output.pickled_trained_model, rs=42),
        threshold = expand(rules.train_test_accuracy_species_prediction_top4_log_reg_thresh.output.threshold, rs=42)
    output:
        predicted_label_df = 'results/tables/{species}_prediction/{species}_predicted_label.tsv',
        scaled_feature_df = 'results/tables/{species}_prediction/{species}_scaled_features.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/predict_species_snoRNA_label_final.py"
