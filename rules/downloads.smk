import os

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
