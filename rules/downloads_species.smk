import os


rule ensembl_species_genome:
    """ Download species genome in FASTA file (.fa) from Ensembl ftp servers."""
    output:
        genome = 'data/references/{species}_genome.fa'
    params:
        link = "ftp://ftp.ensembl.org/pub/release-105/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule ensembl_species_gtf:
    """ Download species genome annotation in gtf file (.gtf) from Ensembl ftp
        servers. While very large, beware that not all organism are present in
        the Ensembl database (examples of species NOT in Ensembl: Anguilla
        anguilla, Xenopus laevis, Drosophila pseudoobscura, Drosophila simulans,
        Branchiostoma lanceolatum, Manis javanica). This mean that, unfortunately,
        you cannot predict snoRNA abundance status with our predictor in these
        species."""
    output:
        gtf = "data/references/{species}.gtf"
    params:
        link = "ftp://ftp.ensembl.org/pub/release-105/gtf/{species}/*[0-9].gtf.gz"
    shell:
        "wget -O temp2.gz {params.link} && "
        "gunzip temp2.gz && "
        "mv temp2 {output.gtf}"


rule rna_central_to_ensembl_id_species:
    """ Download the table to convert RNACentral ids to Ensembl ids and keep
        only snoRNAs of specific species."""
    input:
        gtf = rules.ensembl_species_gtf.output.gtf
    output:
        conversion_table = "data/references/rna_central_to_ensembl_id_{species}.tsv"
    params:
        link = "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings",
        file = lambda wildcards: ensembl_data_dict[wildcards.species]
    shell:
        "read id < <(grep -oE 'gene_id \"([[:upper:]]*[[:lower:]]*|[[:lower:]]*[[:upper:]])' {input.gtf} | head -n1 | sed 's/gene_id \"//') && "
        "wget -O tempfile {params.link}/{params.file} && "
        "grep --color=never $id tempfile | grep --color=never snoRNA > {output.conversion_table} && "
        "rm tempfile"


rule get_taxid_species:
    """ Get the taxid (taxon ID) for a given organism using the tool taxoniq.
        You must give the organism scientific name."""
    output:
        taxid = 'data/references/taxid_{species}.tsv'
    params:
        species_complete_name = lambda wildcards: species_dict[wildcards.species]
    shell:
        """pip3 install 'taxoniq==0.6.0' && """
        """taxoniq --scientific-name '{params.species_complete_name}' url > temp_id && """
        """IN=$(cat temp_id) && arr=(${{IN//id=/ }}) && """
        """echo ${{arr[1]}} | sed s'/\"//' > {output.taxid} && """
        """rm temp_id"""

rule get_species_HG_expression:
    """ Get the expression datasets from normal samples in the Bgee database
        (Bastian et al., NAR, 2021) for the given species. These datasets are
        available for 29 model animal species and processed using a common
        bioinformatic pipeline."""
    output:
        expression_dir = directory("data/references/{species}_Bgee_expression_datasets")
    params:
        species_complete_name = lambda wildcards: species_uppercase_dict[wildcards.species],
        link = "https://bgee.org/ftp/current/download/processed_expr_values/rna_seq"
    shell:
        "wget {params.link}/{params.species_complete_name}/{params.species_complete_name}_RNA-Seq_read_counts_TPM_FPKM.tar.gz && "
        "tar -xf {params.species_complete_name}_RNA-Seq_read_counts_TPM_FPKM.tar.gz && "
        "mkdir -p data/references && "
        "mv {params.species_complete_name}_RNA-Seq_read_counts_TPM_FPKM {output.expression_dir} && "
        "for file in {output.expression_dir}/*; do gunzip $file; done && "
        "rm {params.species_complete_name}_RNA-Seq_read_counts_TPM_FPKM.tar.gz"


rule get_RNA_central_snoRNAs_species:
    """ Get all entries in RNA central marked as snoRNA for a given taxid. Be
        patient, this rule might take up to 30 min depending on your internet
        connection."""
    input:
        taxid = rules.get_taxid_species.output.taxid
    output:
        RNA_central_snoRNAs = 'data/references/rna_central_all_{species}_snoRNAs.tsv'
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
