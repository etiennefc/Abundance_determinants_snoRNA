#!/usr/bin/python3
import pandas as pd
import functions as ft

""" Create a density plot per confusion value comparison (e.g. FP vs TP) for all
    numerical features in the top 10 predictive features per snoRNA type C/D vs
    H/ACA). Each comparison counts only one time a snoRNA (ex: it considers a TP
    snoRNA once even if it is predicted multiple time as a TP across iterations)."""

sno_type = snakemake.wildcards.sno_type
sno_type = sno_type[0] + '/' + sno_type[1:]
color_dict = snakemake.params.color_dict
output = snakemake.output.density
numerical_feature = snakemake.wildcards.top_10_numerical_features
comparison = snakemake.wildcards.comparison_confusion_val
feature_df = pd.read_csv(snakemake.input.feature_df, sep='\t')

# Select only one snoRNA type
feature_df = feature_df[feature_df[sno_type] == 1.0]
feature_df = feature_df[['gene_id_sno', numerical_feature]]


# Get the list of all snoRNAs of a given confusion_value inside dict
sno_per_confusion_value_paths = snakemake.input.sno_per_confusion_value
sno_per_confusion_value = {}
for path in sno_per_confusion_value_paths:
    confusion_value = path.split('/')[-1]
    confusion_value = confusion_value.split('_')[0]
    df = pd.read_csv(path, sep='\t')
    sno_list = df['gene_id_sno'].to_list()
    sno_per_confusion_value[confusion_value] = sno_list


# Get the snoRNA feature value for each confusion value in the comparison
confusion_val1, confusion_val2_3 = comparison.split('_vs_')
confusion_val2, confusion_val3 = confusion_val2_3.split('_')  # this is respectively TN and TP
df1 = feature_df[feature_df['gene_id_sno'].isin(sno_per_confusion_value[confusion_val1])]
df2 = feature_df[feature_df['gene_id_sno'].isin(sno_per_confusion_value[confusion_val2])]
df3 = feature_df[feature_df['gene_id_sno'].isin(sno_per_confusion_value[confusion_val3])]

# Get only snoRNAs that are always predicted as their confusion value
# (i.e. remove snoRNAs that are for example predicted in an iteration as FP and in another as TN)
if confusion_val1 == 'FP':
    all_fp = df1['gene_id_sno'].to_list()
    all_tn = df2['gene_id_sno'].to_list()
    all_tp = df3['gene_id_sno'].to_list()
    all_fn = feature_df[feature_df['gene_id_sno'].isin(sno_per_confusion_value['FN'])]
    all_fn = list(pd.unique(all_fn['gene_id_sno']))
    real_fp = list(set(all_fp) - set(all_tn))
    real_tn = list(set(all_tn) - set(all_fp))
    real_tp = list(set(all_tp) - set(all_fn))
    df1 = df1[df1['gene_id_sno'].isin(real_fp)]
    df2 = df2[df2['gene_id_sno'].isin(real_tn)]
    df3 = df3[df3['gene_id_sno'].isin(real_tp)]
elif confusion_val1 == 'FN':
    all_fn = df1['gene_id_sno'].to_list()
    all_tn = df2['gene_id_sno'].to_list()
    all_tp = df3['gene_id_sno'].to_list()
    all_fp = feature_df[feature_df['gene_id_sno'].isin(sno_per_confusion_value['FP'])]
    all_fp = list(pd.unique(all_fp['gene_id_sno']))
    real_fn = list(set(all_fn) - set(all_tp))
    real_tn = list(set(all_tn) - set(all_fp))
    real_tp = list(set(all_tp) - set(all_fn))
    df1 = df1[df1['gene_id_sno'].isin(real_fn)]
    df2 = df2[df2['gene_id_sno'].isin(real_tn)]
    df3 = df3[df3['gene_id_sno'].isin(real_tp)]

len_df1, len_df2, len_df3 = len(df1), len(df2), len(df3)

# Create density plot
colors = [color_dict[confusion_val1], color_dict[confusion_val2], color_dict[confusion_val3]]
ft.density_x([df1[numerical_feature], df2[numerical_feature], df3[numerical_feature]],
            numerical_feature, 'Density', 'linear', f'{comparison} ({sno_type})',
            colors, [f'{confusion_val1} ({len_df1})', f'{confusion_val2} ({len_df2})', f'{confusion_val3} ({len_df3})'], output)
