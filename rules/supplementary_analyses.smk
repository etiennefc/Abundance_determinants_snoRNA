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
        sno_conservation = config['path']['sno_conservation'],
        intersection_upstream_sno_conservation = config['path']['intersection_upstream_sno_conservation'],
        upstream_sno_conservation = config['path']['upstream_sno_conservation']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/sno_conservation.py"

rule AQR_binding:
    """ Find the overlap between AQR (IBP160) and intron containing snoRNAs """
    input:
        aqr_HepG2_bed = rules.eclip_encode_download_AQR.output.bed_file_1,
        aqr_K562_bed = rules.eclip_encode_download_AQR.output.bed_file_2,
        sno_bed = rules.gtf_to_bed.output.all_sno_bed,
        HG_bed = rules.filter_HG_bed.output.HG_bed,
        df = rules.merge_feature_df.output.feature_df
    output:
        intron_bed = 'data/bed_files/intron_of_all_intronic_snoRNAs.bed',
        overlap_sno_AQR = 'data/references/eCLIP/overlap_sno_AQR.filtered.bed'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/AQR_binding.py"

rule DKC1_binding_haca:
    """ Find the overlap between DKC1 binding (PAR-CLIP and eCLIP) and H/ACA snoRNAs."""
    input:
        dkc1_HepG2_eCLIP_bed = rules.eclip_encode_download_DKC1.output.bed_file,
        dkc1_HEK293_par_clip_bed = rules.par_clip_download.output.bed_dkc1,
        sno_bed = rules.gtf_to_bed.output.all_sno_bed,
        df = rules.merge_feature_df.output.feature_df
    output:
        overlap_sno_DKC1_eCLIP = 'data/references/eCLIP/overlap_sno_DKC1_eCLIP_filtered.bed',
        overlap_sno_DKC1_par_clip = 'data/references/PAR_CLIP/overlap_sno_DKC1_PAR_CLIP_filtered.bed'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/DKC1_binding_eCLIP.py"

rule core_RBP_binding_cd_par_clip:
    """ Find the overlap between NOP58 RepA/RepB, NOP56, FBL/FBL_mnase binding (PAR-CLIP) and C/D snoRNAs."""
    input:
        RBP_par_clip_beds = rules.merge_bed_peaks.output.merge_beds,
        sno_bed = rules.gtf_to_bed.output.all_sno_bed,
        df = rules.merge_feature_df.output.feature_df
    output:
        overlap_sno_NOP58_par_clip = 'data/references/PAR_CLIP/overlap_sno_NOP58_PAR_CLIP_filtered.bed',
        overlap_sno_NOP56_par_clip = 'data/references/PAR_CLIP/overlap_sno_NOP56_PAR_CLIP_filtered.bed',
        overlap_sno_FBL_par_clip = 'data/references/PAR_CLIP/overlap_sno_FBL_PAR_CLIP_filtered.bed'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/core_RBP_binding_cd_par_clip.py"

