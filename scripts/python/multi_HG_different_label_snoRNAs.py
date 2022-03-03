#!/usr/bin/python3
import pandas as pd
import collections as coll

host_df = pd.read_csv(snakemake.input.host_df)
host_df = host_df[['sno_id', 'host_id', 'host_name']]
host_df = host_df.rename(columns={'sno_id': 'gene_id_sno'})
feature_df = pd.read_csv(snakemake.input.feature_df, sep='\t')

# Drop intergenic snoRNAs and merge to host df
feature_df = feature_df[feature_df['abundance_cutoff_host'] != 'intergenic']
feature_df =  feature_df.merge(host_df, how='left', on='gene_id_sno')

same_label, different_label = [], []
same, diff = 0, 0
for i, group in enumerate(feature_df.groupby('host_id')):
    grouped_df = group[1]
    if len(grouped_df) > 1:  # multi-HG
        if len(list(pd.unique(grouped_df['abundance_cutoff_2']))) == 1:  # same label for snoRNAs in same HG
            same_label.append(grouped_df)
            same +=1
        if len(list(pd.unique(grouped_df['abundance_cutoff_2']))) == 2:  # different label for snoRNAs in same HG
            different_label.append(grouped_df)
            diff += 1
print(f'Same label multi-HG: {same}')
print(f'Different label multi-HG: {diff}')

different_label_df = pd.concat(different_label)
different_label_df.to_csv(snakemake.output.multi_HG_different_label_snoRNAs, sep='\t')
same_label_df = pd.concat(same_label)
