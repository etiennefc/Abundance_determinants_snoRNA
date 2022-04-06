import os

""" Quantify snoRNA expression in mouse embryonic stem cells using the
    TGIRT-Seq pipeline."""

include: "downloads.smk"

rule coco_ca:
    """ Generate corrected annotation from the gtf."""
    input:
        gtf = rules.ensembl_mouse_gtf.output.gtf
    output:
        gtf_corrected = config['path']['corrected_mouse_gtf']
    conda:
        "../envs/coco.yaml"
    shell:
        "python3 git_repos/coco/bin/coco.py ca {input.gtf} -o {output.gtf_corrected}"


rule qc_before_trim:
    """Assess fastq quality before trimming reads"""
    input:
        fastq1 = "data/references/mouse_fastq/{id}_1.fastq.gz",
        fastq2 = "data/references/mouse_fastq/{id}_2.fastq.gz"
    output:
        qc_report1 = "data/FastQC/Before_trim/{id}_1_fastqc.html",
        qc_report2 = "data/FastQC/Before_trim/{id}_2_fastqc.html"
    threads:
        32
    params:
        out_dir = "data/FastQC/Before_trim"
    log:
        "log/fastqc/before_trim/{id}.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "-f fastq "
        "-t {threads} "
        "-o {params.out_dir} "
        "{input.fastq1} "
        "{input.fastq2} "
        "&> {log}"


rule trimming:
    """Trims the input FASTQ files using Trimmomatic"""
    input:
        fastq1 = "data/references/mouse_fastq/{id}_1.fastq.gz",
        fastq2 = "data/references/mouse_fastq/{id}_2.fastq.gz"
    output:
        fastq1 = "data/Trimmomatic/trimmed_reads/{id}_R1.fastq.gz",
        fastq2 = "data/Trimmomatic/trimmed_reads/{id}_R2.fastq.gz",
        unpaired_fastq1 = "data/Trimmomatic/trimmed_reads/{id}_R1.unpaired.fastq.gz",
        unpaired_fastq2 = "data/Trimmomatic/trimmed_reads/{id}_R2.unpaired.fastq.gz"
    threads:
        32
    params:
        options = [
            "ILLUMINACLIP:data/Trimmomatic/Adapters-PE_NextSeq.fa:2:12:10:8:true",
            "TRAILING:30", "LEADING:30", "MINLEN:20"
        ]
    log:
        "log/trimmomatic/{id}.log"
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.fastq1} {input.fastq2} "
        "{output.fastq1} {output.unpaired_fastq1} "
        "{output.fastq2} {output.unpaired_fastq2} "
        "{params.options} "
        "&> {log}"


rule qc_after_trim:
    """Assess fastq quality after trimming reads"""
    input:
        fastq1 = rules.trimming.output.fastq1,
        fastq2 = rules.trimming.output.fastq2
    output:
        qc_report1 = "data/FastQC/After_trim/{id}_R1_fastqc.html",
        qc_report2 = "data/FastQC/After_trim/{id}_R2_fastqc.html"
    threads:
        32
    params:
        out_dir = "data/FastQC/After_trim"
    log:
        "log/fastqc/after_trim/{id}.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "-f fastq "
        "-t {threads} "
        "-o {params.out_dir} "
        "{input.fastq1} "
        "{input.fastq2} "
        "&> {log}"


rule star_index:
    """Generate the genome index needed for STAR alignment"""
    input:
        fasta = config['path']['mouse_genome'],
        standard_gtf = config['path']['mouse_gtf']
    output:
        chrNameLength = "data/star_index/chrNameLength.txt"
    threads:
        32
    params:
        index_dir = "data/star_index/"
    log:
        "log/star/star_index.log"
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.standard_gtf} "
        "--sjdbOverhang 74 "
        "&> {log}"


rule star_align:
    """Align reads to reference genome using STAR"""
    input:
        fastq1 = rules.trimming.output.fastq1,
        fastq2 = rules.trimming.output.fastq2,
        idx = rules.star_index.output
    output:
        bam = "results/star/{id}/Aligned.sortedByCoord.out.bam"
    threads:
        32
    params:
        outdir = "results/star/{id}/",
        index_dir = "data/star_index/"
    log:
        "log/star/star_align_{id}.log"
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.index_dir} "
        "--readFilesIn {input.fastq1} {input.fastq2} "
        "--runThreadN {threads} "
        "--readFilesCommand zcat "
        "--outReadsUnmapped Fastx "
        "--outFilterType BySJout "
        "--outStd Log "
        "--outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.outdir} "
        "--outFilterScoreMinOverLread 0.3 "
        "--outFilterMatchNminOverLread 0.3 "
        "--outFilterMultimapNmax 100 "
        "--winAnchorMultimapNmax 100 "
        "--alignEndsProtrude 5 ConcordantPair"
        "&> {log}"


rule coco_cc:
    """Quantify the number of counts, counts per million (CPM) and transcript
        per million (TPM) for each gene using CoCo correct_count (cc)."""
    input:
        gtf = rules.coco_ca.output.gtf_corrected,
        bam = rules.star_align.output.bam
    output:
        counts = "results/coco_cc_mouse/{id}.csv"
    threads:
        32
    params:
        coco_path = "git_repos/coco/bin"
    log:
        "log/coco/coco_{id}.log"
    conda:
        "../envs/coco_cc.yaml"
    shell:
        "python {params.coco_path}/coco.py cc "
        "--countType both "
        "--thread {threads} "
        "--strand 1 "
        "--paired "
        "{input.gtf} "
        "{input.bam} "
        "{output.counts} "
        "&> {log}"

rule merge_coco_output_mouse:
    """ Merge CoCo correct count outputs into one count, cpm or tpm file (all
        samples merged inside one dataframe). This rule takes in input the
        output result directory of CoCo within the TGIRT-Seq pipeline."""
    input:
        counts = expand(rules.coco_cc.output.counts, id=all_sample_ids)
    output:
        merged_counts = os.path.join(config['path']['coco_merge_mouse'], "merged_counts.tsv"),
        merged_cpm = os.path.join(config['path']['coco_merge_mouse'], "merged_cpm.tsv"),
        merged_tpm = os.path.join(config['path']['coco_merge_mouse'], "merged_tpm.tsv")
    params:
        input_dir = ""
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/merge_coco_cc_output_mouse.py"
