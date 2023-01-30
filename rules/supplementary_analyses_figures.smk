rule density_upstream_conservation:
    """ Create density plot of the average conservation of 100 nt upstream of either 
        recent or conserved expressed interenic snoRNAs. """
    input:
        sno_cons = rules.sno_conservation.output.sno_conservation,
        upstream_sno_cons = rules.sno_conservation.output.upstream_sno_conservation,
	feature_label_df = rules.merge_feature_df.output.feature_df
    output:
        density = os.path.join(config['figures']['density'], 'upstream_cons_intergenic_sno.svg')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/density_upstream_conservation.py"

rule stacked_bar_labels_snotype_HG_biotype:
    """ Generate a grouped stacked bar per expressed or not expressed
        snoRNAs of sno_type and host biotype (one stacked bar per hue)."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        bar = os.path.join(config['figures']['bar'],
                            'ab_status_sno_type_host_biotype.svg')
    params:
        host_biotype_colors = config['colors_complex']['host_biotype2'],
        sno_type_colors = config['colors_complex']['sno_type']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/stacked_bar_labels_snotype_HG_biotype.py"

rule dist_to_bp_thresh_other_features:
    """ Look for differences between expressed snoRNAs that are close (<=100nt) vs far (>100nt)
        from the branch point, with regards to terminal stem, box score and target type."""
    input:
        df = rules.merge_feature_df.output.feature_df
    output:
        density_terminal_stem_cd = os.path.join(config['figures']['density'],
                            'dist_to_bp_thresh_terminal_stem_cd.svg'),
        density_terminal_stem_haca = os.path.join(config['figures']['density'],
                            'dist_to_bp_thresh_terminal_stem_haca.svg'),
        density_box_score_cd = os.path.join(config['figures']['density'],
                            'dist_to_bp_thresh_box_score_cd.svg'),
        density_box_score_haca = os.path.join(config['figures']['density'],
                            'dist_to_bp_thresh_box_score_haca.svg'),
        density_sno_mfe_cd = os.path.join(config['figures']['density'],
                            'dist_to_bp_thresh_sno_mfe_cd.svg'),
        density_sno_mfe_haca = os.path.join(config['figures']['density'],
                            'dist_to_bp_thresh_sno_mfe_haca.svg'),
	bar_target_cd = os.path.join(config['figures']['bar'],
                            'dist_to_bp_thresh_sno_target_cd.svg'),
        bar_target_haca = os.path.join(config['figures']['bar'],
                            'dist_to_bp_thresh_sno_target_haca.svg')
    params:
        dist_to_bp_group_colors = config['colors_complex']['dist_to_bp_group'],
        sno_target_colors = config['colors_complex']['sno_target']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/dist_to_bp_thresh_other_features.py"

rule dist_to_bp_thresh_HG_AQR_overlap:
    """ Look for differences between expressed snoRNAs that are close (<=100nt) vs far (>100nt)
        from the branch point, with regards to the binding of AQR (IBP160) in their intron."""
    input:
        df = rules.merge_feature_df.output.feature_df,
        aqr_overlap_HG = rules.AQR_binding.output.overlap_sno_AQR
    output:
        bar_HG_cd = os.path.join(config['figures']['bar'],
                            'dist_to_bp_thresh_HG_AQR_overlap_cd.svg'),
        bar_HG_haca = os.path.join(config['figures']['bar'],
                            'dist_to_bp_thresh_HG_AQR_overlap_haca.svg')
    params:
        AQR_binding_colors = config['colors_complex']['AQR_binding']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/dist_to_bp_thresh_HG_AQR_overlap.py"


rule core_protein_binding_ab_status_bar:
    """ Look for differences between expressed and not expressed snoRNAs
        with regards to the binding of different core proteins (DKC1 for H/ACA; NOP58/56 and FBL for C/D)."""
    input:
        df = rules.merge_feature_df.output.feature_df,
        dkc1_eclip_overlap = rules.DKC1_binding_haca.output.overlap_sno_DKC1_eCLIP,
        dkc1_par_clip_overlap = rules.DKC1_binding_haca.output.overlap_sno_DKC1_par_clip,
        nop58_par_clip_overlap = rules.core_RBP_binding_cd_par_clip.output.overlap_sno_NOP58_par_clip,
        nop56_par_clip_overlap = rules.core_RBP_binding_cd_par_clip.output.overlap_sno_NOP56_par_clip,
        fbl_par_clip_overlap = rules.core_RBP_binding_cd_par_clip.output.overlap_sno_FBL_par_clip
    output:
        bar_haca_dkc1_eclip = os.path.join(config['figures']['bar'],
                            'DKC1_binding_eCLIP_ab_status_haca.svg'),
        bar_haca_dkc1_par_clip = os.path.join(config['figures']['bar'],
                            'DKC1_binding_PAR_CLIP_ab_status_haca.svg'),
        bar_cd_nop58_par_clip = os.path.join(config['figures']['bar'],
                            'NOP58_binding_PAR_CLIP_ab_status_cd.svg'),
        bar_cd_fbl_par_clip = os.path.join(config['figures']['bar'],
                            'FBL_binding_PAR_CLIP_ab_status_cd.svg'),
        bar_cd_nop56_par_clip = os.path.join(config['figures']['bar'],
                            'NOP56_binding_PAR_CLIP_ab_status_cd.svg')
    params:
        RBP_binding_colors = config['colors_complex']['RBP_binding']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/graphs/core_protein_binding_ab_status_bar.py"


