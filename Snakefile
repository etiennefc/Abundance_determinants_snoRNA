import os
from get_figures import get_figures_path

configfile: "config.json"

simple_untreated_id = list(config['mouse_untreated_fastq_ids'].keys())
simple_RA_treated_id = list(config['mouse_RA_treated_fastq_ids'].keys())
all_sample_ids = simple_untreated_id + simple_RA_treated_id
organism_name="Mus musculus"

wildcard_constraints:
    feature_hue = "({})".format("|".join(config["feature_hue"])),
    numerical_features = "({})".format("|".join(config["numerical_features"])),
    numerical_features_scaled = "({})".format("|".join(config["numerical_features_scaled"])),
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
    rs = "({})".format("|".join([str(rs_) for rs_ in config['rs']]))

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
#include: "rules/cv_train_test.smk"
include: "rules/figures_model_output.smk"
include: "rules/other_figures.smk"
#include: "rules/cv_train_test_wo_clusters.smk"
#include: "rules/cv_train_test_wo_feature.smk"
#include: "rules/cv_train_test_one_feature.smk"
#include: "rules/cv_train_test_scale_after_split.smk"
#include: "rules/shap_retest.smk"
include: "rules/cv_train_test_10_iterations.smk"
include: "rules/cv_train_test_10_iterations_only_hamming.smk"
include: "rules/cv_train_test_manual_split.smk"
include: "rules/cv_train_test_manual_split_top3.smk"
include: "rules/cv_train_test_manual_split_top4.smk"
#include: "rules/cv_train_test_20_iterations.smk"
include: "rules/clip_seq.smk"
include: "rules/candidate_analyses.smk"
#include: "rules/kallisto.smk"
include: "rules/shap_subgroups.smk"
include: "rules/tgirt_seq_mouse.smk"
include: "rules/mouse_prediction.smk"
include: "rules/mouse_prediction_figures.smk"
include: "rules/cv_train_test_for_species_prediction_top4.smk"
include: "rules/cv_train_test_for_species_prediction_top4_random_state.smk"
include: "rules/cv_train_test_for_species_prediction_top4_log_reg_thresh.smk"
include: "rules/snora81_overexpression_analyses.smk"



rule all:
    input:
        #### Do 10 different iterations of the cv-train-test
        #test_accuracy_scale_after_split = expand(os.path.join(config['path']['test_accuracy'],
        #                            '{models2}_test_accuracy_scale_after_split_{iteration}.tsv'),
        #                            models2=config['models2'], iteration=config['iteration']),
        #confusion_matrix_f1_scale_after_split = expand(os.path.join(config['path']['confusion_matrix_f1'],
        #                '{models2}_confusion_matrix_w_f1_score_scale_after_split_{iteration}.tsv'),
        #                models2=config['models2'], iteration=config['iteration']),
        # Bed for VAP
        #expressed_intronic_sno_bed = os.path.join(config['path']['sno_bed_vap'],
        #                            "expressed_intronic_snoRNA.bed"),
        # Modify SHAP _beeswarm script to sort by median and not by average or sum of SHAP values
        #fake_log = "log/modify_shap.log",
        #concat_df = config['path']['all_feature_rank_df_10_iterations'],
        #
        # Do 20 different iterations of the cv-train-test
        #test_accuracy_scale_after_split_20 = expand(os.path.join(config['path']['test_accuracy'],
        #                            '{models2}_test_accuracy_scale_after_split_{iteration_20}.tsv'),
        #                            models2=config['models2'], iteration_20=config['iteration_20']),
        #confusion_matrix_f1_scale_after_split_20 = expand(os.path.join(config['path']['confusion_matrix_f1'],
        #                '{models2}_confusion_matrix_w_f1_score_scale_after_split_{iteration_20}.tsv'),
        #                models2=config['models2'], iteration_20=config['iteration_20']),
        #not_expressed_haca = config['path']['not_expressed_haca_fa'],
        #small_expressed_cd = config['path']['small_expressed_cd_fa'],
        #features_df = config['path']['feature_df'],
        #abundance_cutoff_df = config['path']['all_RNA_biotypes_tpm_df_cutoff'],
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

        # Train with only combined_box_hamming, sno_mfe and terminal_stem_mfe as features
        training_accuracy_top3 = expand(os.path.join(config['path']['training_accuracy'],
                                            '{models2}_training_accuracy_scale_after_split_top3_{manual_iteration}.tsv'), **config),
        test_accuracy_top3 = expand(os.path.join(config['path']['test_accuracy'],
                                            '{models2}_test_accuracy_scale_after_split_top3_{manual_iteration}.tsv'), **config),

        # Train with only combined_box_hamming, sno_mfe, terminal_stem_mfe and host_expressed as features
        training_accuracy_top4 = expand(os.path.join(config['path']['training_accuracy'],
                                            '{models2}_training_accuracy_scale_after_split_top4_{manual_iteration}.tsv'), **config),
        test_accuracy_top4 = expand(os.path.join(config['path']['test_accuracy'],
                                            '{models2}_test_accuracy_scale_after_split_top4_{manual_iteration}.tsv'), **config),
        confusion_matrix_top4 = expand(os.path.join(config['path']['confusion_matrix_f1'],
                                            '{models2}_confusion_matrix_w_f1_score_scale_after_split_top4_{manual_iteration}.tsv'), **config),

        #sno_presence_test_sets = config['path']['sno_presence_test_sets'],
        #sno_per_confusion_value = expand(os.path.join(config['path']['sno_per_confusion_value'], '{confusion_value}_snoRNAs_10_iterations.tsv'), **config),
        #concat_shap = expand(os.path.join(config['path']['shap_10_iterations'], '{models2}_{confusion_value}_shap_values_concat.tsv'), **config),
        #merge_beds = expand(os.path.join(config['path']['par_clip'], '{rbp}_merged.bed'),
        #                    rbp=['DKC1', 'FBL', 'FBL_mnase', 'NOP56', 'NOP58_repA', 'NOP58_repB']),
        #real_confusion_value_df = expand(os.path.join(config['path']['real_confusion_value'], '{confusion_value}_w_features.tsv'), **config),
        #transcript_tpm = "results/kallisto/transcript_tpm.tsv",
        mapped_snoRNA_bed = expand(os.path.join(config['path']['par_clip'], '{rbp}_mapped_to_snoRNAs.bed'),
                            rbp=['NOP58_repA', 'NOP58_repB', 'NOP56', 'FBL', 'FBL_mnase', 'DKC1']),
        multi_HG_different_label_snoRNAs = config['path']['multi_HG_different_label_snoRNAs'],
        sno_per_confusion_value = config['path']['sno_per_confusion_value_manual_split'],


        # Mouse snoRNA quantification
        qc_before_trim = expand("data/FastQC/Before_trim/{id}_1_fastqc.html", id=all_sample_ids),
        qc_after_trim = expand("data/FastQC/After_trim/{id}_R1_fastqc.html", id=all_sample_ids),
        #coco_cc_mouse = os.path.join(config['path']['coco_merge_mouse'], "merged_tpm.tsv"),

        feature_df_mouse = config['path']['feature_df_mouse'],
        test_accuracy_mouse = expand(os.path.join(config['path']['test_accuracy_mouse'],'{models2}_test_accuracy_{manual_iteration}.tsv'), **config),
        confusion_value_mouse = os.path.join(config['path']['confusion_matrix_f1_mouse'], 'consensus_confusion_value_per_sno.tsv'),
        consensus_conf_val_df_per_model = expand(os.path.join(config['path']['confusion_matrix_f1_mouse'], 'consensus_confusion_value_per_sno_{models2}.tsv'), **config),

        # Species prediction
        test_accuracy_species_prediction = expand(os.path.join(config['path']['test_accuracy_mouse'],
                                                '{models2}_test_accuracy_top4_species_prediction.tsv'), **config),
        confusion_matrix_species_prediction = expand(os.path.join(config['path']['confusion_matrix_f1_mouse'],
                                    '{models2}_confusion_matrix_w_f1_score_top4_species_prediction.tsv'), **config),

        # Species prediction 5 random_state datasets split
        training_accuracy_rs = expand(os.path.join(config['path']['training_accuracy'],
                                    '{models2}_training_accuracy_top4_species_prediction_{rs}.tsv'), **config),
        #test_accuracy_rs = expand(os.path.join(config['path']['test_accuracy_mouse'],
        #                            '{models2}_test_accuracy_top4_species_prediction_{rs}.tsv'), **config),

        # Log reg thresh species prediction
        confusion_matrix_log_reg_thresh = expand(os.path.join(config['path']['confusion_matrix_f1_mouse'],
                        'log_reg_confusion_matrix_w_f1_score_top4_species_prediction_thresh_{rs}.tsv'), **config),

        # SNORA81 overexpression analyses
        #predicted_label = expand('results/tables/snora81_overexpression/{models2}_SNORA81_label_{manual_iteration}.tsv', **config)

        #fake_log2 = 'log/fake_log_snora77b.log'
        #hamming_distance_box_df = config['path']['hamming_distance_box_df'],
        #h_aca_box_location_expressed = config['path']['h_aca_box_location_expressed'],
        #a = config['path']['c_d_and_prime_box_location_not_expressed'],
        ##log_create_env = "log/create_local_env.log",  #to put in all_downloads
        ##log_create_env_2 = "log/create_local_env_2.log",
        #one_hot_encoded_df = config['path']['one_hot_encoded_df_before_split'],
        ##HG_not_expressed_sno_bed = os.path.join(config['path']['sno_bed_vap'],
        ##                            "HG_not_expressed_snoRNA.bed"),
        ##test_accuracy = expand(os.path.join(config['path']['test_accuracy'], ##Added knn with models2
        ##                '{models2}_test_accuracy.tsv'), models2=config['models2']),
        ##confusion_matrix = expand(os.path.join(config['path']['confusion_matrix_f1'],
        ##                '{models}_confusion_matrix_w_f1_score.tsv'), models=config['models']),
        ##top_features_df = config['path']['all_feature_rank_df'],
        ##shap_local_FN_log_scale_after_split = expand("log/shap_local_FN_{models2}_scale_after_split.log", models2=config['models2']),
        ##shap_local_FP_log_scale_after_split = expand("log/shap_local_FP_{models2}_scale_after_split.log", models2=config['models2']),
        ##shap_local_TN_log_scale_after_split = expand("log/shap_local_TN_{models2}_scale_after_split.log", models2=config['models2']),
        ##shap_local_TP_log_scale_after_split = expand("log/shap_local_TP_{models2}_scale_after_split.log", models2=config['models2']),
        ##snora77b_terminal_stem = config['path']['snora77b_terminal_stem_fasta'],
        ##test_accuracy_wo_clusters = expand(os.path.join(config['path']['test_accuracy'],
        ##                            '{models2}_test_accuracy_wo_clusters.tsv'), models2=config['models2']),
        ##confusion_matrix_wo_clusters = expand(os.path.join(config['path']['confusion_matrix_f1'],
        ##                '{models}_confusion_matrix_w_f1_score_wo_clusters.tsv'), models=config['models']),
        ##test_accuracy_wo_feature_effect = expand(os.path.join(config['path']['test_accuracy'],
        ##                            '{models2}_test_accuracy_wo_{feature_effect}.tsv'),
        ##                            models2=config['models2'], feature_effect=config['feature_effect']),
        ##confusion_matrix_wo_conservation = expand(os.path.join(config['path']['confusion_matrix_f1'],
        ##                '{models}_confusion_matrix_w_f1_score_wo_{feature_effect}.tsv'),
        ##                models=config['models'], feature_effect=config['feature_effect']),
        ##test_accuracy_one_feature = expand(os.path.join(config['path']['test_accuracy'],
        ##                            '{models2}_test_accuracy_only_{one_feature}.tsv'),
        ##                            models2=config['models2'], one_feature=config['one_feature']),
        ##confusion_matrix_one_feature = expand(os.path.join(config['path']['confusion_matrix_f1'],
        ##                '{models}_confusion_matrix_w_f1_score_only_{one_feature}.tsv'),
        ##                models=config['models'], one_feature=config['one_feature']),

        ### Do scaling after split to avoid overestimate of the models performance
        ##test_accuracy_scale_after_split = expand(os.path.join(config['path']['test_accuracy'],
        ##                            '{models2}_test_accuracy_scale_after_split.tsv'), models2=config['models2']),
        ##confusion_matrix_f1_scale_after_split = expand(os.path.join(config['path']['confusion_matrix_f1'],
        ##                '{models2}_confusion_matrix_w_f1_score_scale_after_split.tsv'), models2=config['models']),

        # Consensus sequence of SNORA77B-terminal stem across vertebrates
        #cs = "CS.stk",


        #bp_distance = config['path']['bp_distance']
        #sno_HG_df = config['path']['sno_tpm_df'],
        #sno_bed_path = config['path']['all_sno_bed']
        #snhg14_sno_bed = os.path.join(config['path']['sno_bed'], "snhg14_snoRNA.bed")

        #snhg14_bed = config['path']['snhg14_bed'],
        #intronic_sno_bed = os.path.join(config['path']['sno_bed'], "all_intronic_snoRNA.bed"),
        #output_table = config['path']['sno_location_table'],
        #hg_gtf = config['path']['hg_gtf'],
        #temp = 'data/temp.txt',
        ##bp_distance = config['path']['bp_distance'],
        #HG_bed = config['path']['hg_bed']
        #sno_bp_distance = config['path']['sno_distance_bp'],
        #sno_conservation = config['path']['sno_conservation'],
        #sno_length = config['path']['sno_length'],
        #mfe_final = config['path']['structure_stability_tsv'],
        #snodb_formatted = config['path']['snodb_formatted'],
        #mfe_stem = config['path']['terminal_stem_mfe_tsv'],
        #length_stem = config['path']['terminal_stem_length'],
        #tpm_cutoff = config['path']['sno_tpm_df_cutoff'],


        ##variable_sno_labels_df = config['path']['variable_sno_label_same_host']

rule all_downloads:
    input:
        snhg14_bed = config['path']['snhg14_bed'],
        lnctard = config['path']['lnctard'],
        eclip_bed = config['path']['TARDBP_rep1_eclip_bed'],
        bed_dkc1 = config['path']['DKC1_par_clip'],
        mouse_untreated_fastq_1_gz = expand("data/references/mouse_fastq/{id_untreated}_1.fastq.gz", id_untreated=simple_untreated_id),
        mouse_RA_treated_fastq_1_gz = expand("data/references/mouse_fastq/{id_RA_treated}_1.fastq.gz", id_RA_treated=simple_RA_treated_id),
        genome = config['path']['mouse_genome'],
        gtf = config['path']['mouse_gtf'],
        coco_git = 'git_repos/coco',  # don't forget to create the git env required for this download
        conversion_table = config['path']['rna_central_to_ensembl_id'],
        RNA_central_snoRNAs = 'data/references/rna_central_all_mouse_snoRNAs.tsv',
        recounts = 'data/recount_datasets.csv'

rule all_figures:
    input:
        # Modify SHAP _beeswarm script to sort by median and not by average or sum of SHAP values
        fake_log = "log/modify_shap.log",
        figures = get_figures_path(config)
