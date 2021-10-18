import os

include: "downloads.smk"
include: "data_processing.smk"


rule bw_to_bg:
    """ Convert phastCons bigwig file into bedgraph file and sort that bedgraph
        in place by chromosome and start. Also remove unwanted tabs (\t) at the
        end of some lines in the bed file of all snoRNAs, add 'chr' at the
        beginning of all lines and sort by chromosome and start."""
    input:
        phastcons_bigwig = rules.phastcons_download.output.phastcons
    output:
        phastcons_bedgraph = config['path']['phastcons_bg']
    params:
        sno_bed_path = config['path']['all_sno_bed']
    conda:
        "../envs/conservation.yaml"
    shell:
        """bigWigToBedGraph {input.phastcons_bigwig} {output.phastcons_bedgraph} && """
        """sort -k1,1 -k2,2n -o {output.phastcons_bedgraph} {output.phastcons_bedgraph} && """
        """sed -i 's/\t$//g; s/^/chr/g' {params.sno_bed_path} && """
        """sort -k1,1 -k2,2n -o {params.sno_bed_path} {params.sno_bed_path}"""

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
        sno_conservation = config['path']['sno_conservation']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/python/sno_conservation.py"
