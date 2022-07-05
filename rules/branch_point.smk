import os

include: 'data_processing.smk'
include: 'downloads.smk'


rule find_bedtools_dir:
    """ Find the exact location of binary bedtools that is needed as an
        input for branchpointer. It is CRUCIAL that only one branchpointer
        R env is created otherwise this step will not work with mutiple R
        envs."""
    params:
        dir = config['path']['r_env_dir']
    output:
        bedtools_dir = "data/bedtools_dir.txt"
    shell:
        "find .snakemake/conda/*/bin/ -name R > {output.bedtools_dir} && "
        "sed -i 's/R/bedtools/g' {output.bedtools_dir}"

rule branch_point:
    """ Calculate the distance between intronic snoRNAs and predicted
        branch point position."""
    input:
        transcript_id_df = rules.sno_exon_location.output.output_table,
        hg_gtf = rules.generate_HG_gtf.output.hg_gtf,
        bedtools_dir = rules.find_bedtools_dir.output.bedtools_dir,
        genome = rules.ensembl_genome.output.genome  
    output:
        bp_window_table = config['path']['bp_window_table'],
        bp_distance = config['path']['bp_distance']
    params:
        bedtools_dir = config['path']['bedtools_dir']
    conda:
        "../envs/branchpointer.yaml"
    script:
        "../scripts/r/branch_pointer.R"

rule get_best_bp:
    """ From all the nucleotides proposed to be the branch point (bp)
        location by the rule branch_point, select only the one with the
        highest predicted probability and return a simplified version of
        the df returned by branch_point. Return also an output table
        containing snoRNA distance to the best branch point per HG
        intron."""
    input:
        bp_distance_total_df = rules.branch_point.output.bp_distance,
        sno_location_df = rules.sno_exon_location.output.output_table
    output:
        bp_distance_simple = config['path']['bp_distance_simple'],
        sno_distance_bp = config['path']['sno_distance_bp']
    params:
        sno_overlap_df = os.path.join(
                            config['path']['sno_location_exon'],
                            "sno_overlap_hg.tsv")
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/get_best_bp.py"
