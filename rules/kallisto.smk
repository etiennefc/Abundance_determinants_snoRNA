import os

include: "downloads.smk"

rule create_transcriptome:
    """ Uses gffread to generate a transcriptome """
    input:
        genome = config['path']['genome_v101'],
        gtf = config['path']['gtf']
    output:
        seqs = config['path']['transcriptome']
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output.seqs}"

rule kallisto_index:
    """ Generates the transcriptome index for Kallisto """
    input:
        transcriptome = rules.create_transcriptome.output.seqs
    output:
        idx = config['path']['kallisto_index']
    params:
        kmer = "31"
    log:
        "log/kallisto/index.log"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index "
        "--index={output.idx} "
        "--kmer-size={params.kmer} "
        "{input.transcriptome} "
        "&> {log}"

rule kallisto_quant:
    """ Generates counts using Kallisto pseudo-alignment """
    input:
        idx = rules.kallisto_index.output.idx,
        fq1 = os.path.join(config['path']['trimmed_reads'], "{tissue}_R1.fastq.gz"),
        fq2 = os.path.join(config['path']['trimmed_reads'], "{tissue}_R2.fastq.gz")
    output:
        quant = "results/kallisto/{tissue}/abundance.tsv",
        h5 = "results/kallisto/{tissue}/abundance.h5"
    params:
        bootstrap = "50",
        outdir = "results/kallisto/{tissue}"
    log:
        "log/kallisto/{tissue}.log"
    threads:
        8
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}"

rule generate_transcriptID_geneName:
    """
    Generating a two-column text file containing the gene -> transcript
    relationship
    """
    input:
        gtf = config['path']['gtf']
    output:
        map = config['path']['gene_name']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/generate_transcriptID_geneName.py"

rule combine_gene_quantification:
    """
    Custom Python script to collect and format Kallisto results for further
    processing.
    """
    input:
        datasets = expand(rules.kallisto_quant.output.quant, **config),
        map = rules.generate_transcriptID_geneName.output.map
    output:
        tpm = "results/kallisto/tpm.tsv",
        est_counts = "results/kallisto/est_counts.tsv",
        transcript_tpm = "results/kallisto/transcript_tpm.tsv",
        transcript_est_counts = "results/kallisto/transcript_est_counts.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/combine_gene_quantification.py"
