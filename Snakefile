import os
from get_figures import get_figures_path

configfile: "config.json"

simple_untreated_id = list(config['mouse_untreated_fastq_ids'].keys())
simple_RA_treated_id = list(config['mouse_RA_treated_fastq_ids'].keys())
all_sample_ids = simple_untreated_id + simple_RA_treated_id
organism_name = "Mus musculus"
species_original_name = ['Rattus norvegicus',  'Danio rerio',
                         'Gallus gallus', 'Ornithorhynchus anatinus',
                        'Xenopus tropicalis', 'Pan troglodytes', 'Bos taurus',
                        'Oryctolagus cuniculus', 'Gorilla gorilla', 'Macaca mulatta']
species = [name.replace(' ', '_').lower() for name in species_original_name]
species_uppercase = [name.replace(' ', '_') for name in species_original_name]
species_dict, species_uppercase_dict, ensembl_data_dict = {}, {}, {}
for i, name in enumerate(species_original_name):
    species_dict[species[i]] = name
    species_uppercase_dict[species[i]] = species_uppercase[i]
    if name  == 'Drosophila melanogaster':
        ensembl_data_dict[species[i]] = 'flybase.tsv'
    elif name  == 'Caenorhabditis elegans':
        ensembl_data_dict[species[i]] = 'wormbase.tsv'
    else:
        ensembl_data_dict[species[i]] = 'ensembl.tsv'


wildcard_constraints:
    feature_hue = "({})".format("|".join(config["feature_hue"])),
    numerical_features = "({})".format("|".join(config["numerical_features"])),
    numerical_features_scaled = "({})".format("|".join(config["numerical_features_scaled"])),
    mouse_numerical_features = "({})".format("|".join(config["mouse_numerical_features"])),
    categorical_features = "({})".format("|".join(config["categorical_features"])),
    intronic_features = "({})".format("|".join(config["intronic_features"])),
    sno_type = "({})".format("|".join(config["sno_type"])),
    hg_biotype = "({})".format("|".join(config["hg_biotype"])),
    FN = "({})".format("|".join(config["FN"])),
    feature_effect = "({})".format("|".join(config["feature_effect"])),
    pca_hue = "({})".format("|".join(config["pca_hue"])),
    models2 = "({})".format("|".join(config["models2"])),
    iteration = "({})".format("|".join(config["iteration"])),
    iteration_20 = "({})".format("|".join(config["iteration_20"])),
    manual_iteration = "({})".format("|".join(config["manual_iteration"])),
    top_10_numerical_features = "({})".format("|".join(config["top_10_numerical_features"])),
    top_10_categorical_features = "({})".format("|".join(config["top_10_categorical_features"])),
    comparison_confusion_val = "({})".format("|".join(config["comparison_confusion_val"])),
    interesting_sno_ids = "({})".format("|".join(config["interesting_sno_ids"])),
    id_untreated = "({})".format("|".join(list(config['mouse_untreated_fastq_ids'].keys()))),
    id_RA_treated = "({})".format("|".join(list(config['mouse_RA_treated_fastq_ids'].keys()))),
    rs = "({})".format("|".join([str(rs_) for rs_ in config['rs']])),
    species = "({})".format("|".join(species))

#Include data processing rules to generate the dataset
include: "rules/data_processing.smk"
include: "rules/downloads.smk"
include: "rules/branch_point.smk"
include: "rules/conservation.smk"
include: "rules/structure.smk"
include: "rules/terminal_stem.smk"
include: "rules/chip_seq.smk"
include: "rules/merge_features.smk"
include: "rules/figures_features.smk"
include: "rules/feature_normalization.smk"
include: "rules/figures_model_output.smk"
include: "rules/other_figures.smk"
include: "rules/cv_train_test_10_iterations.smk"
include: "rules/cv_train_test_10_iterations_only_hamming.smk"
include: "rules/cv_train_test_manual_split.smk"
include: "rules/cv_train_test_manual_split_top3.smk"
include: "rules/cv_train_test_manual_split_top4.smk"
include: "rules/clip_seq.smk"
include: "rules/candidate_analyses.smk"
include: "rules/shap_subgroups.smk"
include: "rules/tgirt_seq_mouse.smk"
include: "rules/mouse_prediction.smk"
include: "rules/mouse_prediction_figures.smk"
include: "rules/yeast_prediction.smk"
include: "rules/cv_train_test_for_species_prediction_top4.smk"
include: "rules/cv_train_test_for_species_prediction_top4_random_state.smk"
include: "rules/cv_train_test_for_species_prediction_top4_log_reg_thresh.smk"
include: "rules/cv_train_test_for_species_prediction_top3_random_state.smk"
include: "rules/cv_train_test_for_species_prediction_top3_log_reg_thresh.smk"
include: "rules/gtex_HG_cutoff.smk"
include: "rules/cv_train_test_manual_split_gtex_HG.smk"
include: "rules/downloads_species.smk"
include: "rules/species_prediction.smk"
include: "rules/species_prediction_figures.smk"


rule all:
    input:
        # Do manual split of 10 test sets
        confusion_matrix = expand(os.path.join(config['path']['confusion_matrix_f1'],
                                '{models2}_confusion_matrix_w_f1_score_scale_after_split_{manual_iteration}.tsv'), **config),
        test_accuracy = expand(os.path.join(config['path']['test_accuracy'],
                                    '{models2}_test_accuracy_scale_after_split_{manual_iteration}.tsv'), **config),
        get_all_shap = expand(os.path.join(config['path']['shap_10_iterations'], '{models2}_{manual_iteration}_shap_values.tsv'), **config),

        # Train with only combined_hamming as feature
        training_accuracy_only_hamming = expand(os.path.join(config['path']['training_accuracy'],
                                            '{models2}_training_accuracy_scale_after_split_only_hamming_{manual_iteration}.tsv'), **config),
        test_accuracy_only_hamming = expand(os.path.join(config['path']['test_accuracy'],
                                            '{models2}_test_accuracy_scale_after_split_only_hamming_{manual_iteration}.tsv'), **config),
        confusion_matrix_only_hamming = expand(os.path.join(config['path']['confusion_matrix_f1'],
                                            '{models2}_confusion_matrix_w_f1_score_scale_after_split_only_hamming_{manual_iteration}.tsv'), **config),
        # Train with only combined_box_hamming, sno_mfe and terminal_stem_mfe as features
        training_accuracy_top3 = expand(os.path.join(config['path']['training_accuracy'],
                                            '{models2}_training_accuracy_scale_after_split_top3_{manual_iteration}.tsv'), **config),
        test_accuracy_top3 = expand(os.path.join(config['path']['test_accuracy'],
                                            '{models2}_test_accuracy_scale_after_split_top3_{manual_iteration}.tsv'), **config),
        confusion_matrix_top3 = expand(os.path.join(config['path']['confusion_matrix_f1'],
                                            '{models2}_confusion_matrix_w_f1_score_scale_after_split_top3_{manual_iteration}.tsv'), **config),
        # Train with only combined_box_hamming, sno_mfe, terminal_stem_mfe and host_expressed as features
        training_accuracy_top4 = expand(os.path.join(config['path']['training_accuracy'],
                                            '{models2}_training_accuracy_scale_after_split_top4_{manual_iteration}.tsv'), **config),
        test_accuracy_top4 = expand(os.path.join(config['path']['test_accuracy'],
                                            '{models2}_test_accuracy_scale_after_split_top4_{manual_iteration}.tsv'), **config),
        confusion_matrix_top4 = expand(os.path.join(config['path']['confusion_matrix_f1'],
                                            '{models2}_confusion_matrix_w_f1_score_scale_after_split_top4_{manual_iteration}.tsv'), **config),
        multi_HG_different_label_snoRNAs = config['path']['multi_HG_different_label_snoRNAs'],
        sno_per_confusion_value_manual_split = config['path']['sno_per_confusion_value_manual_split']
	

rule all_mouse:                                                                                                                                                                                                              input:
	# Mouse snoRNA quantification
        #qc_before_trim = expand("data/FastQC/Before_trim/{id}_1_fastqc.html", id=all_sample_ids),
        #qc_after_trim = expand("data/FastQC/After_trim/{id}_R1_fastqc.html", id=all_sample_ids),
        coco_cc_mouse = os.path.join(config['path']['coco_merge_mouse'], "merged_tpm.tsv"),
        feature_df_mouse = config['path']['feature_df_mouse'],
        test_accuracy_mouse = expand(os.path.join(config['path']['test_accuracy_mouse'],'{models2}_test_accuracy_{manual_iteration}.tsv'), **config),
        confusion_value_mouse = os.path.join(config['path']['confusion_matrix_f1_mouse'], 'consensus_confusion_value_per_sno.tsv'),
        consensus_conf_val_df_per_model = expand(os.path.join(config['path']['confusion_matrix_f1_mouse'], 'consensus_confusion_value_per_sno_{models2}.tsv'), **config),



rule all_downloads:
    input:
        snhg14_bed = config['path']['snhg14_bed'],
        tpm_df = os.path.join(config['path']['merge_coco_output'], 'tpm_v101.csv'),
        human_gtf = config['path']['gtf'],
        human_gtf_df = config['path']['gtf_tsv_table'],
        human_genome_fa = config['path']['genome_v101'],
        snodb = config['path']['snodb'],
        lnctard = config['path']['lnctard'],
        hg_df = config['path']['host_gene_df'],
        mouse_untreated_fastq_1_gz = expand("data/references/mouse_fastq/{id_untreated}_1.fastq.gz", id_untreated=simple_untreated_id),
        mouse_RA_treated_fastq_1_gz = expand("data/references/mouse_fastq/{id_RA_treated}_1.fastq.gz", id_RA_treated=simple_RA_treated_id),
        mouse_genome = config['path']['mouse_genome'],
        mouse_gtf = config['path']['mouse_gtf'],
        coco_git = 'git_repos/coco',  # don't forget to create the git env required for this download
        conversion_table = config['path']['rna_central_to_ensembl_id'],
        RNA_central_snoRNAs = 'data/references/rna_central_all_mouse_snoRNAs.tsv',
        recounts = 'data/recount_datasets.csv',
        gtex_data = 'data/references/gtex_tpm_df.tsv',
        gtex_data_unpaired = 'data/references/gtex_unpaired_tpm_df.tsv'

rule all_figures:
    input:
        # Modify SHAP _beeswarm script to sort by median and not by average or sum of SHAP values
        fake_log = "log/modify_shap.log",
        figures = get_figures_path(config)


rule species_downloads:
    input:
        genome = expand('data/references/{species}_genome.fa', species=species),
        gtf = expand("data/references/{species}.gtf", species=species),
        conversion_table = expand("data/references/rna_central_to_ensembl_id_{species}.tsv", species=species),
        RNA_central_species_snoRNAs = expand('data/references/rna_central_all_{species}_snoRNAs.tsv', species=species),
        expression_dir = expand("data/references/{species}_Bgee_expression_datasets", species=species)

rule species_predictions:
    input:
        predicted_snoRNA_labels = expand('results/tables/{species}_prediction/{species}_predicted_label.tsv', species=species)

rule species_figures:
    input:
        bar_ab_status_prediction_species = os.path.join(config['figures']['bar'],
                                                'ab_status_prediction_all_species.svg'),
        scatter_species_prediction_sno_nb = os.path.join(config['figures']['scatter'],
                    'sno_nb_ab_status_prediction_all_species.svg'),
        df = 'results/tables/summary_table_sno_type_host_biotype_species.tsv'
