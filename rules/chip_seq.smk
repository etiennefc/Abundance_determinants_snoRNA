import os

include: "data_processing.smk"


rule bed_for_vap:
    """ Generate bed files that will be used by the tool VAP for bedgraph
        aggregation. One bed file for expressed intronic snoRNAs, one bed for
        not expressed intronic snoRNAs, one bed for expressed intergenic snoRNAs,
        one bed for not expressed intergenic snoRNAs, one bed for HG of expressed
        snoRNAs and one bed for HG of not expressed snoRNAs."""
    input:
        intronic_sno_bed = rules.generate_snoRNA_beds.output.intronic_sno_bed,
        intergenic_sno_bed = rules.generate_snoRNA_beds.output.intergenic_sno_bed,
        hg_bed = rules.HG_gtf_to_bed.output.HG_bed,
        hg_df = config['path']['host_gene_df'],
        abundance_status_df = config['path']['feature_df']
    output:
        expressed_intronic_sno_bed = os.path.join(config['path']['sno_bed_vap'],
                                    "expressed_intronic_snoRNA.bed"),
        not_expressed_intronic_sno_bed = os.path.join(config['path']['sno_bed_vap'],
                                    "not_expressed_intronic_snoRNA.bed"),
        expressed_intergenic_sno_bed = os.path.join(config['path']['sno_bed_vap'],
                                    "expressed_intergenic_sno.bed"),
        not_expressed_intergenic_sno_bed = os.path.join(config['path']['sno_bed_vap'],
                                    "not_expressed_intergenic_sno.bed"),
        HG_expressed_sno_bed = os.path.join(config['path']['sno_bed_vap'],
                        "HG_expressed_snoRNA.bed"),
        HG_not_expressed_sno_bed = os.path.join(config['path']['sno_bed_vap'],
                                    "HG_not_expressed_snoRNA.bed")
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/bed_for_vap.py"


#rule vap_aggregate_beds:

#rule extract_avg_peaks_sno:
#    """ Extract average value in the ~ 500 pb upstream of all snoRNAs."""

#rule extract_avg_peaks_HG:
#    """ Extract average value in the ~ 500/1000 pb upstream of the first exon of
#        intronic snoRNA HGs."""
