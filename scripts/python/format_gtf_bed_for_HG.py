#!/usr/bin/python3
import pandas as pd
import subprocess as sp

df = pd.read_csv(snakemake.input.gtf_bed, sep='\t',
                    names=['chr', 'start', 'end', 'gene_id', 'dot',
                            'strand', 'source', 'feature', 'dot2', 'attributes'], index_col=False)

embedded_genes = ['miRNA', 'Mt_tRNA', 'scaRNA', 'scRNA', 'snoRNA', 'snRNA', 'sRNA']
embedded_genes = [f'gene_biotype "{item}"' for item in embedded_genes]
embedded_genes = '{}'.format('|'.join(embedded_genes))

# Keep only gene features and remove embedded genes
df = df[df['feature'] == 'gene']
df = df[~df['attributes'].str.contains(embedded_genes)]

# Add "chr" in front of chr number
df['chr'] = 'chr' + df['chr'].astype(str)
df.to_csv('temp_gtf.bed', index=False, sep='\t', header=False)

# Sort gtf bed file
sp.call('sort -k1,1 -k2,2n temp_gtf.bed > '+snakemake.output.formatted_gtf_bed, shell=True)
sp.call('rm temp_gtf.bed', shell=True)
