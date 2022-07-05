import os

rule coco_human_tpm_df:
    """ Download TGIRT-Seq abundance dataset of all human genes quantified by CoCo."""
    output:
        tpm_df = os.path.join(
                    config['path']['merge_coco_output'],
                    'tpm_v101.csv')
    params:
        link = config['download']['coco_human_tpm_df']
    shell:
        "wget -O {output.tpm_df} {params.link}"

rule human_gtf_download:
    """ Download the annotation (gtf file) of all human genes."""
    output:
        gtf = config['path']['gtf']
    params:
        link = config['download']['gtf_human']
    shell:
        "wget -O {output.gtf} {params.link}"

rule gtf_tsv_table_download:
    """ Download the tsv table (gtf converted to tsv) of all human genes."""
    output:
        gtf_df = config['path']['gtf_tsv_table']
    params:
        link = config['download']['gtf_tsv_table']
    shell:
        "wget -O {output.gtf_df} {params.link}"

rule refseq_gtf:
    """ Download Refseq annotation from ftp server and generate a gtf of only the
        SNHG14 gene. Generate also a bed file out of this SNHG14 gtf. The "15"
        in the following lines refer to chr15 in which SNHG14 is encoded"""
    output:
        snhg14_gtf = config['path']['snhg14_gtf'],
        snhg14_bed = config['path']['snhg14_bed']
    params:
        link = config['download']['refseq_gtf']
    shell:
        """wget {params.link} -O refseq_temp.gtf.gz --quiet && zcat refseq_temp.gtf.gz | """
        """grep SNHG14 | sed 's/NC_000015.10/15/g' > {output.snhg14_gtf} && """
        """zcat refseq_temp.gtf.gz | grep SNHG14 | """
        """awk -v OFS='\t' 'NR>1 {{print "chr15", $4, $5, "ENSG00000224078", """
        """$6, $7, $2, $3, $8, $1, $12, "TEST"$0"TEST2", $10, $16, "transcript_name", $16}}' | """
        """sed -E 's/;//g; s/TEST.*exon_number..//g; s/..TEST2//g' > {output.snhg14_bed} && """
        """rm refseq_temp.gtf.gz"""

rule ensembl_genome:
    """ Download human genome in FASTA file (.fa) from Ensembl ftp servers and
        remove scaffolds (starting with KI270728.1) at the end of the file."""
    output:
        genome = config['path']['genome_v101']
    params:
        link = config['download']['ensembl_genome']
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "sed '/>KI270728.1/Q' temp > {output.genome} && "
        "rm temp"

rule phastcons_download:
    """ Download phastCons conservation score across 100 vertebrates for all
        nucleotides in the human genome (hg38) in a bigwig format."""
    output:
        phastcons = config['path']['phastcons_bw']
    params:
        link = config['download']['phastcons']
    shell:
        "wget -O {output.phastcons} {params.link}"

rule snodb_nmd_di_promoters_download:
    """ Download the snodb table used for this analysis, and the NMD substrates
        and dual-initiation promoter host genes tables from Zenodo."""
    output:
        snodb = config['path']['snodb'],
        nmd = config['path']['nmd_lykke_andersen_corrigendum'],
        di_promoter = config['path']['di_promoter_nepal']
    params:
        link_snodb = config['download']['snodb'],
        link_nmd = config['download']['nmd_lykke_andersen_corrigendum'],
        link_di_promoter = config['download']['di_promoter_nepal']
    shell:
        "wget -O {output.snodb} {params.link_snodb} && "
        "wget -O {output.snodb} {params.link_nmd} && "
        "wget -O {output.snodb} {params.link_di_promoter}"

rule host_gene_list_download:
    """ Download the host genes (and functions) of human snoRNAs from 
        reference table in Zenodo."""
    output:
        hg_df = config['path']['host_gene_df'],
        hg_pc_function = config['path']['protein_coding_HG_functions']
    params:
        link_hg = config['download']['final_host_gene_list_v101'],
        link_pc_functions = config['download']['protein_coding_HG_functions']
    shell:
        "wget -O {output.hg_df} {params.link_hg} && "
        "wget -O {output.hg_pc_function} {params.link_pc_functions}"

rule lncTarD_download:
    """ Download the table of lncRNA (non-coding host genes) validated functions
        in human diseases."""
    output:
        lnctard = config['path']['lnctard']
    params:
        link_lnctard = config['download']['lnctard']
    shell:
        "wget -O temp.txt {params.link_lnctard} && "
        "cut -f 21 temp.txt | sort -r | uniq > {output.lnctard} && "
        "rm temp.txt"

rule forgi_viennrna_download:
    """ Download the folding tool called forgi from the ViennRNA package."""
    output:
        forgi_log = "log/forgi.log"
    shell:
        "pip install --upgrade forgi &> {output.forgi_log}"

rule eclip_encode_download:
    """ Download the eCLIP-Seq datasets for the TARDBP RNA binding protein.
        Convert the bigBed into a bed file."""
    output:
        bed_file_1 = config['path']['TARDBP_rep1_eclip_bed'],
        bed_file_2 = config['path']['TARDBP_rep2_eclip_bed']
    params:
        link_bed_file_1 = config['download']['TARDBP_rep1_eclip'],
        link_bed_file_2 = config['download']['TARDBP_rep2_eclip']
    conda:
        "../envs/bigbedtobed.yaml"
    shell:
        "wget -O temp_bigbed_1.bb {params.link_bed_file_1} && "
        "bigBedToBed temp_bigbed_1.bb {output.bed_file_1} && "
        "wget -O temp_bigbed_2.bb {params.link_bed_file_2} && "
        "bigBedToBed temp_bigbed_2.bb {output.bed_file_2} && "
        "rm temp_bigbed_*"

rule liftover_chain_file_download:
    """ Download the chain file (hg19 to hg38) needed for the liftover of par-clip beds"""
    output:
        liftover_chain_file = config['path']['liftover_chain_file']
    params:
        link_chain_file = config['download']['liftover_chain_file']
    shell:
        "wget -O {output.liftover_chain_file}.gz {params.link_chain_file} && "
        "gunzip {output.liftover_chain_file}"

rule par_clip_download:
    """ Download the PAR-CLIP datasets of NOP58, NOP56, FBL and DKC1 from the
        Kishore et al. 2013 Genome Biology paper. Sort them by chr and start."""
    input:
        chain_file = rules.liftover_chain_file_download.output.liftover_chain_file
    output:
        bed_nop58_repA = config['path']['NOP58_repA_par_clip'],
        bed_nop58_repB = config['path']['NOP58_repB_par_clip'],
        bed_nop56 = config['path']['NOP56_par_clip'],
        bed_fbl = config['path']['FBL_par_clip'],
        bed_fbl_mnase = config['path']['FBL_mnase_par_clip'],
        bed_dkc1 = config['path']['DKC1_par_clip']
    params:
        link_bed_nop58_repA = config['download']['NOP58_repA_par_clip'],
        link_bed_nop58_repB = config['download']['NOP58_repB_par_clip'],
        link_bed_nop56 = config['download']['NOP56_par_clip'],
        link_bed_fbl = config['download']['FBL_par_clip'],
        link_bed_fbl_mnase = config['download']['FBL_mnase_par_clip'],
        link_bed_dkc1 = config['download']['DKC1_par_clip']
    conda:
        "../envs/liftover.yaml"
    shell:
        "paths=$(echo {params}) && "
        "arr=(${{paths// / }}) && "
        "outputs=$(echo {output}) && "
        "arr_outputs=(${{outputs// / }}) && "
        "for index in ${{!arr[@]}}; do "
        "echo $index; "
        "echo ${{arr[$index]}}; "
        "wget -O temp_par${{index}}.gz ${{arr[$index]}}; "
        "gunzip temp_par${{index}}.gz; "
        "liftOver temp_par${{index}} {input.chain_file} temp_liftover${{index}} unmapped_${{index}}; "
        "sort -k1,1 -k2,2n temp_liftover${{index}} > ${{arr_outputs[$index]}}; "
        "done; "
        "rm temp_par* && rm temp_liftover* && rm unmapped_*"

rule sashimi_script_download:
    """ Download the script needed to create sashimi plots."""
    output:
        sashimi_script = config['path']['sashimi_script']
    params:
        link = config['download']['sashimi_script']
    shell:
        "wget {params.link} -O {output.sashimi_script} && chmod u+x {output.sashimi_script}"

rule mouse_tgirt_seq_untreated_fastq_downloads:
    """ Download the 6 fastq from the small RNA-Seq (with TGIRT) experiment done
        in McCann et al. 2020 (NAR). These are mouse untreated embryonic stem
        cells (n=3)."""
    output:
        mouse_fastq_1_gz = "data/references/mouse_fastq/{id_untreated}_1.fastq.gz",
        mouse_fastq_2_gz = "data/references/mouse_fastq/{id_untreated}_2.fastq.gz",
    params:
        link = lambda wildcards: config['mouse_untreated_fastq_ids'][wildcards.id_untreated]
    shell:
        "wget {params.link}.1.fastq.gz.1 -O {output.mouse_fastq_1_gz} && "
        "wget {params.link}.2.fastq.gz.1 -O {output.mouse_fastq_2_gz}"

rule mouse_tgirt_seq_RA_treated_fastq_downloads:
    """ Download the 6 fastq from the small RNA-Seq (with TGIRT) experiment done
        in McCann et al. 2020 (NAR). These are mouse embryonic stem cells
        treated with retinoic acid (RA) (n=3)."""
    output:
        mouse_fastq_1_gz = "data/references/mouse_fastq/{id_RA_treated}_1.fastq.gz",
        mouse_fastq_2_gz = "data/references/mouse_fastq/{id_RA_treated}_2.fastq.gz",
    params:
        link = lambda wildcards: config['mouse_RA_treated_fastq_ids'][wildcards.id_RA_treated]
    shell:
        "wget {params.link}.1.fastq.gz.1 -O {output.mouse_fastq_1_gz} && "
        "wget {params.link}.2.fastq.gz.1 -O {output.mouse_fastq_2_gz}"

rule ensembl_mouse_genome:
    """ Download mouse genome in FASTA file (.fa) from Ensembl ftp servers and
        remove scaffolds (starting with JH584299.1) at the end of the file."""
    output:
        genome = config['path']['mouse_genome']
    params:
        link = config['download']['ensembl_mouse_genome']
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "sed '/>JH584299.1/,$d; s/>/>chr/' temp > {output.genome} && "
        "rm temp"

rule ensembl_mouse_gtf:
    """ Download mouse genome annotation in gtf file (.gtf) from Ensembl ftp
        servers."""
    output:
        gtf = config['path']['mouse_gtf']
    params:
        link = config['download']['ensembl_mouse_gtf']
    shell:
        "wget -O temp2.gz {params.link} && "
        "gunzip temp2.gz && "
        "mv temp2 {output.gtf}"

rule install_pairedBamToBed12:
    output:
        directory("scripts/pairedBamToBed12/")
    params:
        pairedBamToBed12_bin = 'scripts/pairedBamToBed12/bin'
    conda:
        "../envs/coco.yaml"
    shell:
        'cd scripts && pwd && '
        'git clone https://github.com/Population-Transcriptomics/pairedBamToBed12 && '
        'cd pairedBamToBed12 && '
        'make '

rule download_coco_git:
    """ Download git repository of CoCo."""
    input:
        paired_bam_to_bed12_dependency = rules.install_pairedBamToBed12.output
    output:
        git_coco_folder = directory('git_repos/coco')
    params:
        git_coco_link = config['download']['coco_git_link']
    conda:
        '../envs/git.yaml'
    shell:
        'mkdir -p {output.git_coco_folder} '
        '&& git clone {params.git_coco_link} {output.git_coco_folder}'

rule rna_central_to_ensembl_id:
    """ Download the table to convert RNACentral ids to Ensembl ids."""
    output:
        conversion_table = config['path']['rna_central_to_ensembl_id']
    params:
        link = config['download']['rna_central_to_ensembl_id']
    shell:
        "wget -O tempfile {params.link} && "
        "grep --color=never ENSMUSG tempfile | grep --color=never snoRNA > {output.conversion_table} && "
        "rm tempfile"

rule get_taxid:
    """ Get the taxid (taxon ID) for a given organism using the tool taxoniq.
        You must give the organism scientific name."""
    output:
        taxid = 'data/references/taxid.tsv'
    shell:
        """pip3 install 'taxoniq==0.6.0' && """
        """taxoniq --scientific-name '{organism_name}' url > temp_id && """
        """IN=$(cat temp_id) && arr=(${{IN//id=/ }}) && """
        """echo ${{arr[1]}} | sed s'/\"//' > {output.taxid} && """
        """rm temp_id"""

rule get_taxid_yeast:
    """ Get the taxid (taxon ID) for S. cerevisiae using the tool taxoniq.
        You must give the organism scientific name."""
    output:
        taxid = 'data/references/taxid_yeast.tsv'
    shell:
        """pip3 install 'taxoniq==0.6.0' && """
        """taxoniq --scientific-name 'Saccharomyces cerevisiae' url > temp_id_sacch && """
        """IN=$(cat temp_id_sacch) && arr=(${{IN//id=/ }}) && """
        """echo ${{arr[1]}} | sed s'/\"//' > {output.taxid} && """
        """rm temp_id_sacch"""

rule dowload_mouse_HG_RNA_seq_datasets:
    """ Download RNA-Seq datasets of mouse tissues (including mESC) using recount3"""
    output:
        dataset = 'data/recount_datasets.csv'
    conda:
        "../envs/recount3.yaml"
    script:
        "../scripts/r/download_mouse_HG_RNA_seq_datasets.R"

rule download_yeast_genome:
    """Download the S. cerevisiae reference genome (fasta file) used for this analysis from
        ENSEMBL ftp servers."""
    output:
        genome = config['path']['yeast_genome']
    params:
        link = config['download']['ensembl_yeast_genome']
    shell:
        "wget --quiet -O temp_yeast_genome.fa.gz {params.link} && "
        "gunzip temp_yeast_genome.fa.gz && "
        "sed 's/>/>chr/' temp_yeast_genome.fa > {output.genome} && "
        "rm temp_yeast_genome.fa"

rule download_yeast_annotation:
    """Download the S. cerevisiae annotation (gtf file) used for this analysis."""
    output:
        std_gtf = config['path']['yeast_gtf']
    params:
        link_std_annotation = config['download']['ensembl_yeast_gtf']
    shell:
        "wget --quiet -O temp_yeast.gtf.gz {params.link_std_annotation} && "
        "gunzip temp_yeast.gtf.gz && mv temp_yeast.gtf {output.std_gtf}"

rule download_gtex_host_data:
    """ Download TPM table from GTEx to compare the host gene abundance to that
        computed in our TGIRT-Seq datasets (paired triplicates for the seven
        tissues i.e. brain, liver, skletal_muscle, ovary, breast, testis and
        prostate). We define the gtex_var only for this bash script."""
    output:
        gtex_data = 'data/references/gtex_tpm_df.tsv'
    params:
        link = config['download']['gtex_data'],
        ids = config['gtex_id']
    shell:
        "gtex_var=$(echo {params.ids}) ./scripts/bash/gtex_download.sh {params.link} {output.gtex_data}"

rule download_gtex_host_data_unpaired:
    """ Download TPM table from GTEx to compare the host gene abundance to that
        computed in our TGIRT-Seq datasets. Use this time 7 unpaired tissues (i.e.
        not present in our TGIRT-Seq datasets: adrenal gland, colon, spleen,
        heart, kidney, thyroid and nerve). We define the gtex_var only for this
        bash script."""
    output:
        gtex_data = 'data/references/gtex_unpaired_tpm_df.tsv'
    params:
        link = config['download']['gtex_data'],
        ids = config['gtex_id_unpaired']
    shell:
        "gtex_var=$(echo {params.ids}) ./scripts/bash/gtex_download.sh {params.link} {output.gtex_data}"

rule get_RNA_central_snoRNAs:
    """ Get all entries in RNA central marked as snoRNA for a given taxid"""
    input:
        taxid = rules.get_taxid.output.taxid
    output:
        RNA_central_snoRNAs = 'data/references/rna_central_all_mouse_snoRNAs.tsv'
    conda:
        "../envs/psql.yaml"
    shell:
        """
        taxid_var=$(cat {input.taxid}) &&
        psql postgres://reader:NWDMCE5xdipIjRrp@hh-pgsql-public.ebi.ac.uk:5432/pfmegrnargs \
        -c "\copy
         (
            SELECT DISTINCT ON (upi, region_start, region_stop, database)
                r.upi,
                p.short_description,
                p.taxid, a.database,
                a.external_id,
                a.optional_id,
                a.gene,
                a.gene_synonym,
                a.description,
                r.len,
                r.seq_short,
                rfam.rfam_model_id
            FROM rna r
            LEFT JOIN
                rnc_rna_precomputed p ON p.upi = r.upi
            LEFT JOIN rnc_sequence_regions s
                ON s.urs_taxid = p.id
            LEFT JOIN xref x
                ON x.upi = r.upi
            LEFT JOIN rnc_accessions a
                ON a.accession = x.ac
            LEFT JOIN rfam_model_hits rfam
                ON rfam.upi = r.upi
            WHERE p.taxid = $taxid_var
                AND a.ncrna_class IN ('snoRNA', 'scaRNA')
         )  TO '{output.RNA_central_snoRNAs}' WITH (FORMAT CSV, DELIMITER E'\t', HEADER)"
        """
