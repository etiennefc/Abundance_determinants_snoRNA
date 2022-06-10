#!/usr/bin/python3
import pandas as pd
import collections as coll
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder

host_biotype_dfs, sno_type_dfs = snakemake.input.host_biotype_df, snakemake.input.snoRNA_type_df

# Load and merge sno_type df, host_biotype df and predicted ab_status for all species
dfs = []
for i, temp_df in enumerate(snakemake.input.df):
    species_name = temp_df.split('/')[-1].split('_predicted_label')[0].replace('_', ' ').capitalize()
    df = pd.read_csv(temp_df, sep='\t')
    df = df[['gene_id_sno', 'predicted_label']]
    host_biotype_df = pd.read_csv(host_biotype_dfs[i], sep='\t')
    df = df.merge(host_biotype_df[['gene_id_sno', 'host_biotype']], how='left', on='gene_id_sno')
    snoRNA_type_df = pd.read_csv(sno_type_dfs[i], sep='\t')
    snoRNA_type_df = snoRNA_type_df.rename(columns={'gene_id': 'gene_id_sno'})
    df = df.merge(snoRNA_type_df[['gene_id_sno', 'snoRNA_type']], how='left', on='gene_id_sno')
    df['species_name'] = species_name
    dfs.append(df)

# Concat vertically all species dfs
concat_df = pd.concat(dfs)

# Create a simplified version of host_biotype column
simplified_biotype_dict = {'lncRNA': 'non_coding', 'protein_coding': 'protein_coding',
                            'TEC': 'non_coding', 'unitary_pseudogene': 'non_coding',
                            'unprocessed_pseudogene': 'non_coding', 'pseudogene': 'non_coding',
                            'processed_pseudogene': 'non_coding', 'polymorphic_pseudogene': 'non_coding',
                            'processed_transcript': 'non_coding', 'lincRNA': 'non_coding',
                            'antisense': 'non_coding', 'sense_intronic': 'non_coding',
                            'sense_overlapping': 'non_coding'}
concat_df['host_biotype2'] = concat_df['host_biotype'].map(simplified_biotype_dict)
concat_df['host_biotype2'] = concat_df['host_biotype2'].fillna('intergenic')

# One hot encode snoRNA_type and host_biotype2 cols
def split_cols(df, col):
    # Split a column in n one-hot-encoded columns where n is the number of
    # different possibilities in said column
    possibilities = list(pd.unique(df[col]))
    for poss in possibilities:
        df.loc[df[col] == poss, poss] = 1
        df[poss] = df[poss].fillna(0)
    df = df.drop(columns=col)
    return df

sno_type = split_cols(concat_df, 'snoRNA_type')
sno_type = sno_type.drop(columns=['host_biotype2', 'host_biotype'])
host = split_cols(concat_df, 'host_biotype2')

# Merge final_df
final_df = sno_type.merge(host[['gene_id_sno', 'intergenic', 'protein_coding', 'non_coding']], how='left', on='gene_id_sno')
final_df = final_df[['species_name', 'predicted_label', 'H/ACA', 'C/D', 'protein_coding', 'non_coding', 'intergenic']]

# Get number and % of given category (1 (one-hot-encoded)) in col of df
def get_percent(df):
    cols = df.columns
    col_temp, percent_temp = [], []
    for col in cols:
        if col == 'predicted_label':
            label = f'{pd.unique(df[col])[0]}'
            col_temp.append(col)
            percent_temp.append(label)
        elif (col != 'species_name') & (col != 'predicted_label'):
            d = dict(coll.Counter(df[col]))
            if 1 not in d.keys():
                d[1] = 0
            nb = d[1]
            percent = round((nb/len(df)) * 100, 1)
            nb_percent = f'{nb} ({percent}%)'
            col_temp.append(col)
            percent_temp.append(nb_percent)
    temp_df = pd.DataFrame([percent_temp], columns=col_temp)
    return temp_df

# Groupby species_name and get the summary (%) of expressed/not_expressed snoRNAs per characteristics
grouped_df = final_df.groupby(['species_name'])
dfs_final = []
for i, group in grouped_df:
    name, total_sno_nb = i, len(group)
    expressed = group[group['predicted_label'] == 'expressed']
    not_expressed = group[group['predicted_label'] == 'not_expressed']
    temp_dfs = []
    for df_ in [expressed, not_expressed]:
        ddf = get_percent(df_)
        ddf['total'] = len(df_)
        temp_dfs.append(ddf)
    concat_temp_df = pd.concat(temp_dfs)
    concat_temp_df['Species name'] = name
    dfs_final.append(concat_temp_df)


merged_final_df = pd.concat(dfs_final)
merged_final_df = merged_final_df[['Species name', 'predicted_label', 'C/D', 'H/ACA', 'protein_coding', 'non_coding', 'intergenic', 'total']]

merged_final_df.to_csv(snakemake.output.df, sep='\t', index=False)
