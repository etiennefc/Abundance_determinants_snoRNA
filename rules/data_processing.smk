import os

include: "downloads.smk"

'''
rule merge_coco_output:
    """ Merge CoCo correct count outputs into one count, cpm or tpm file (all
        tissues merged inside one dataframe). This rule takes in input the
        output result directory of CoCo within the TGIRT-Seq pipeline."""
    input:
        input_dir = config['path']['coco_output_dir']
    output:
        output_dir = config['path']['merge_coco_output']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/merge_coco_cc_output.py"
'''

rule add_gene_biotype:
    """ Add gene biotype column to tpm file using reference table."""
    input:
        tpm_df = os.path.join(
                    config['path']['merge_coco_output'],
                    'tpm_v101.csv'),
        ref_table = config['path']['gtf_tsv_table']
    output:
        tpm_biotype = os.path.join(
                    config['path']['merge_coco_output'],
                    'tpm_w_biotype_v101.csv')
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/add_gene_biotype.py"

rule add_HG_tpm:
    """ Add host gene abundance columns (for each tissue) for each snoRNA within
        a snoRNA dataframe."""
    input:
        tpm_biotype_df = rules.add_gene_biotype.output.tpm_biotype,
        ref_HG_table = config['path']['host_gene_df']
    output:
        sno_HG_df = config['path']['sno_tpm_df']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/add_HG_tpm.py"

rule abundance_cutoff:
    """ Define a cutoff to classify snoRNAs as 'expressed' or 'not_expressed'
        across tissues (this is the label that will be used by the predictor).
        Define the same cutoff for the host genes (used as another feature:
        HG expressed, not expressed or no HG). """
    input:
        tpm_df = rules.add_HG_tpm.output.sno_HG_df
    output:
        abundance_cutoff_df = config['path']['sno_tpm_df_cutoff']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/abundance_cutoff.py"

rule abundance_cutoff_all_biotypes:
    """ Define a cutoff to classify all types of RNA as 'expressed' or
        'not_expressed' across tissues."""
    input:
        tpm_df = rules.add_gene_biotype.output.tpm_biotype
    output:
        abundance_cutoff_df = config['path']['all_RNA_biotypes_tpm_df_cutoff']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/abundance_cutoff_all_biotypes.py"    

rule generate_HG_gtf:
    """ Generate a gtf containing only the host genes (including SNHG14) from a
        complete gtf. The SNHG14 gtf is generate via the rule refseq_gtf."""
    input:
        host_info_df = config['path']['host_gene_df'],
        gtf = config['path']['gtf'],
        snhg14_gtf = rules.refseq_gtf.output.snhg14_gtf
    output:
        hg_gtf = config['path']['hg_gtf']
    shell:
        """awk -v FS="," 'NR> 1 {{print $8}}' {input.host_info_df} | sort | uniq > hg_temp && """
        """grep -Ff hg_temp {input.gtf} | grep -v SNHG14 > hg_temp_gtf && """
        """cat hg_temp_gtf {input.snhg14_gtf} > {output.hg_gtf} && """
        """rm hg_temp && rm hg_temp_gtf"""

rule gtf_to_bed:
    """ Convert a gtf file into a bed file (first command) and generate out of
        it a bed file of all snoRNA genes (second command)."""
    input:
        gtf = config['path']['gtf']
    output:
        gtf_bed = config['path']['sorted_gtf_bed'],
        all_sno_bed = config['path']['all_sno_bed']
    shell:
        """awk -v OFS='\t' 'NR>6 {{print $1, $4, $5, "to_remove"$10"to_remove", $6, $7, $2, $3, $8, "to_delete"$0}}' {input.gtf} | """
        """sed -E 's/to_remove"//g; s/";to_remove//g; s/to_delete.*gene_id/gene_id/g' | """
        """sort -n -k1,1 -k2,2 > {output.gtf_bed} && """
        """awk '$8=="gene" {{print $0}}' {output.gtf_bed}  | grep snoRNA | sed 's/\t$//g; s/^/chr/g' | sort -k1,1 -k2,2n > {output.all_sno_bed}"""

rule HG_gtf_to_bed:
    """ Convert a HG gtf into a bed format."""
    input:
        HG_gtf = rules.generate_HG_gtf.output.hg_gtf
    output:
        HG_bed = config['path']['hg_bed']
    shell:
        """awk -v OFS='\t' 'NR>6 {{print $1, $4, $5, "to_remove"$10"to_remove", $6, $7, $2, $3, $8, "to_delete"$0}}' {input.HG_gtf} | """
        """sed -E 's/to_remove"//g; s/";to_remove//g; s/to_delete.*gene_id/gene_id/g' | """
        """sort -n -k1,1 -k2,2 > {output.HG_bed} """

rule generate_snoRNA_beds:
    """ From a bed file containing all snoRNAs (generated via gtf_to_bed.sh) and
        a csv file containing snoRNA info (i.e. their HG), generate snoRNA bed
        files (intronic, intergenic, intronic without SNHG14 snoRNAs, only
        SNHG14 snoRNAs)"""
    input:
        all_sno_bed = rules.gtf_to_bed.output.all_sno_bed,
        sno_info_df = rules.add_HG_tpm.output.sno_HG_df
    output:
        intronic_sno_bed = os.path.join(config['path']['sno_bed'], "all_intronic_snoRNA.bed"),
        intergenic_sno_bed = os.path.join(config['path']['sno_bed'], "all_intergenic_snoRNA.bed"),
        intronic_sno_bed_wo_snhg14 = os.path.join(config['path']['sno_bed'], "all_intronic_snoRNA_wo_snhg14.bed"),
        snhg14_sno_bed = os.path.join(config['path']['sno_bed'], "snhg14_snoRNA.bed")
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/generate_snoRNA_beds.py"

rule sno_exon_location:
    """ Get the intron number and distance to closest exons of each intronic
        snoRNA."""
    input:
        sorted_gtf_bed = rules.HG_gtf_to_bed.output.HG_bed,
        sno_tpm_df = rules.add_HG_tpm.output.sno_HG_df,
        sno_HG_coordinates = config['path']['host_gene_df'],
        sno_bed_wo_snhg14 = rules.generate_snoRNA_beds.output.intronic_sno_bed_wo_snhg14,
        snhg14_bed = rules.refseq_gtf.output.snhg14_bed,
        sno_snhg14_bed = rules.generate_snoRNA_beds.output.snhg14_sno_bed
    params:
        sno_location_exon_dir = config['path']['sno_location_exon']
    output:
        output_table = config['path']['sno_location_table']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/sno_location_exon.py"

rule sno_length:
    """ Get the snoRNA length for all snoRNAs (intronic and intergenic)."""
    input:
        all_sno_bed = rules.gtf_to_bed.output.all_sno_bed
    output:
        sno_length = config['path']['sno_length']
    shell:
        """awk -v OFS="\t" '{{print $4, $3-$2+1}}' {input.all_sno_bed} > {output.sno_length}"""

rule HG_functions:
    """ Associate protein-coding and non-coding host genes with their relative
        function. For protein-coding HGs, the functions were manually curated
        from Uniprot; for non-coding HGs (lncRNAs), the functions were retrieved
        from lncTarD."""
    input:
        nc_functions = rules.lncTarD_download.output.lnctard,
        host_genes = config['path']['host_gene_df'],
        pc_functions = config['path']['protein_coding_HG_functions']
    output:
        hg_function_df = config['path']['host_functions']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/host_functions.py"


rule format_snodb:
    """ Remove duplicates from snodb table and format column names. Add host
        gene info from 3 tables (hg biotype and function, NMD susceptibility and
        presence of dual-initiation promoters)."""
    input:
        snodb_original_table = rules.snodb_nmd_di_promoters_download.output.snodb,
        host_gene_df = config['path']['host_gene_df'],
        host_functions = rules.HG_functions.output.hg_function_df,
        nmd = rules.snodb_nmd_di_promoters_download.output.nmd,
        di_promoter = rules.snodb_nmd_di_promoters_download.output.di_promoter
    output:
        snodb_formatted = config['path']['snodb_formatted']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/snodb_table_formatting.py"
