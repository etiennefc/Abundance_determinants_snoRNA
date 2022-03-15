#!/usr/bin/python3
import pandas as pd
import functions as ft
import collections as coll

sno_per_confusion_value_df = pd.read_csv(snakemake.input.sno_per_confusion_value, sep='\t')
host_df = pd.read_csv(snakemake.input.host_df)
multi_HG_different_label_snoRNAs_df = pd.read_csv(snakemake.input.multi_HG_different_label_snoRNAs, sep='\t')
feature_df = pd.read_csv(snakemake.input.feature_df, sep='\t')
feature_df = feature_df[['gene_id_sno', 'abundance_cutoff_2']]
color_dict = snakemake.params.color_dict
bar_output = snakemake.output.hbar


# Drop intergenic snoRNAs and merge sno_per_confusion_value_df to host_df
sno_per_confusion_value_df = sno_per_confusion_value_df.merge(feature_df, how='left', on='gene_id_sno')
intronic_sno_df = host_df.merge(sno_per_confusion_value_df, how='left', left_on='sno_id', right_on='gene_id_sno')



# Select only multi_HG (drop HG with only 1 snoRNA)
mono_HG, multi_HG = [], []
for i, group in enumerate(intronic_sno_df.groupby('host_id')):
    grouped_df = group[1]
    host_id = group[0]
    if len(grouped_df) == 1:  # mono-HG
        mono_HG.append(host_id)
    elif len(grouped_df) > 1:  # multi-HG
        multi_HG.append(host_id)

intronic_sno_df = intronic_sno_df[intronic_sno_df['host_id'].isin(multi_HG)]


# Find in the multi_HG those that have snoRNAs with all the same label vs those with different labels
sno_ids_multi_HG_diff_labels = list(multi_HG_different_label_snoRNAs_df.gene_id_sno)
intronic_sno_df.loc[~intronic_sno_df['gene_id_sno'].isin(sno_ids_multi_HG_diff_labels), 'label_type'] = 'same_label'
intronic_sno_df.loc[intronic_sno_df['gene_id_sno'].isin(sno_ids_multi_HG_diff_labels), 'label_type'] = 'diff_label'

# Find for multi_HG with same snoRNA labels the % of all_expressed or all_not_expressed snoRNAs
all_expressed, all_not_expressed = [], []
for i, group in enumerate(intronic_sno_df[intronic_sno_df['label_type'] == 'same_label'].groupby('host_id')):
    grouped_df = group[1]
    host_id = group[0]
    if 'expressed' in list(grouped_df.abundance_cutoff_2):
        all_expressed.append(host_id)
    elif 'not_expressed' in list(grouped_df.abundance_cutoff_2):
        all_not_expressed.append(host_id)
intronic_sno_df.loc[intronic_sno_df['host_id'].isin(all_expressed), 'expression_category'] = 'all_sno_expressed'
intronic_sno_df.loc[intronic_sno_df['host_id'].isin(all_not_expressed), 'expression_category'] = 'all_sno_not_expressed'


# Find for multi_HG with different snoRNA labels the % of snoRNAs that are 50-50 expressed-not_expressed, more expressed, or more not_expressed
half_expressed_not_expressed, more_expressed, more_not_expressed = [], [], []
for i, group in enumerate(intronic_sno_df[intronic_sno_df['label_type'] == 'diff_label'].groupby('host_id')):
    grouped_df = group[1]
    host_id = group[0]
    if len(grouped_df[grouped_df['abundance_cutoff_2'] == 'expressed']) == len(grouped_df[grouped_df['abundance_cutoff_2'] == 'not_expressed']):
        half_expressed_not_expressed.append(host_id)
    elif len(grouped_df[grouped_df['abundance_cutoff_2'] == 'expressed']) > len(grouped_df[grouped_df['abundance_cutoff_2'] == 'not_expressed']):
        more_expressed.append(host_id)
    elif len(grouped_df[grouped_df['abundance_cutoff_2'] == 'expressed']) < len(grouped_df[grouped_df['abundance_cutoff_2'] == 'not_expressed']):
        more_not_expressed.append(host_id)

intronic_sno_df.loc[intronic_sno_df['host_id'].isin(half_expressed_not_expressed), 'expression_category'] = 'half_expressed_not_expressed'
intronic_sno_df.loc[intronic_sno_df['host_id'].isin(more_expressed), 'expression_category'] = 'more_expressed'
intronic_sno_df.loc[intronic_sno_df['host_id'].isin(more_not_expressed), 'expression_category'] = 'more_not_expressed'



# Groupby label_type, expression_category, host_id and then number of sno per HG (in descending order)
groupby_df = intronic_sno_df.groupby(['label_type', 'expression_category', 'host_name']).size().reset_index()
groupby_df.columns = ['label_type', 'expression_category', 'host_name', 'number_of_sno_per_HG']
groupby_df = groupby_df.sort_values(['label_type', 'expression_category', 'number_of_sno_per_HG'], ascending=[True, True, False])
original_index = list(groupby_df.host_name)
inverted_index = original_index[::-1]  # we invert the order so that it appears correctly on the hbar chart (from top to bottom instead of botton to top)

# Count the number of sno per confusion value per HG
counts = []
for host_name in inverted_index:
    temp = []
    for confusion_val in ['TN', 'TP', 'FN', 'FP']:
        df = intronic_sno_df[intronic_sno_df['host_name'] == host_name]
        number = len(df[df['confusion_matrix'] == confusion_val])
        temp.append(number)
    counts.append(temp)

# Create df to generate the hbar chart (keep the good index determined by the groupby)
final_df = pd.DataFrame(counts, index=inverted_index, columns=['TN', 'TP', 'FN', 'FP'])
print(final_df)
ft.barh(final_df, 10, 20, '', 'Number of snoRNAs', 'Host gene name', bar_output, stacked=True, color=color_dict, width=0.85)
