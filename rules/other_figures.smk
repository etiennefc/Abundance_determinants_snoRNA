import os

rule bar_abundance_status_per_biotype:
    """ Create a stacked bar chart of the percentage of experssed or
        not_expressed RNAs per main RNA biotype (snoRNA, snRNA, tRNA, lncRNA,
        protein-coding)."""
    input:
        abundance_cutoff_df = rules.abundance_cutoff_all_biotypes.output.abundance_cutoff_df
    output:
        bar = os.path.join(config['figures']['bar'],
                    'abundance_status_per_biotype.svg')
    conda:
        "../envs/python.yaml"
    params:
        colors = config['colors_complex']['label']
    script:
        "../scripts/python/graphs/bar_abundance_status_per_biotype.py"
