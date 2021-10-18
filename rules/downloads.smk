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
