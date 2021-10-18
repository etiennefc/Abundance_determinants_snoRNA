#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from upsetplot import plot as upset

""" Generate an upset plot per confusion value (true positive, false positive,
    false negative or true negative) to see the intersection in the snoRNAs
    (mis)classified by all four models."""
df_output_path = "results/tables/confusion_matrix_f1/merged_confusion_matrix.tsv"
vals = ["TP", "TN", "FP", "FN"]
val_cols = ['confusion_matrix_val_log_reg', 'confusion_matrix_val_svc',
            'confusion_matrix_val_gbm', 'confusion_matrix_val_knn']
vals.remove(snakemake.wildcards.confusion_value)

log_reg = pd.read_csv(snakemake.input.log_reg, sep='\t')
log_reg = log_reg[['gene_id_sno', 'confusion_matrix_val_log_reg']]
svc = pd.read_csv(snakemake.input.svc, sep='\t')
svc = svc[['gene_id_sno', 'confusion_matrix_val_svc']]
gbm = pd.read_csv(snakemake.input.gbm, sep='\t')
gbm = gbm[['gene_id_sno', 'confusion_matrix_val_gbm']]
knn = pd.read_csv(snakemake.input.knn, sep='\t')
knn = knn[['gene_id_sno', 'confusion_matrix_val_knn']]

# Merge all dfs and create one df per confusion_value wildcard
df = log_reg.merge(svc, how='left', left_on='gene_id_sno', right_on='gene_id_sno')
df = df.merge(gbm, how='left', left_on='gene_id_sno', right_on='gene_id_sno')
df = df.merge(knn, how='left', left_on='gene_id_sno', right_on='gene_id_sno')
df.to_csv(df_output_path, sep='\t', index=False)


# Select rows containing at least one confusion_value wildcard (ex: TP)
val_df = df[(df.iloc[:, 1:4] == snakemake.wildcards.confusion_value).any(axis=1)]

# Convert confusion_value (ex: TP) to True and all other possible wildcards values (ex: TN, FN, FP) to False
val_df = val_df.replace([snakemake.wildcards.confusion_value, "({})".format("|".join(vals))],
                        [True, False], regex=True)

# Generate Multi_index Serie from DataFrame
upset_df = val_df.copy()
upset_df['gene_id_sno'] = 1  # This ensures that we can sum over the gene_id_sno column with the following groupby
upset_df = upset_df.set_index(val_cols)  # Create multi index

a = upset_df.groupby(level=val_cols).sum()  # Groupby multi index
a = a['gene_id_sno']  # Convert the one-column dataframe into a series

# Create the upset plot
plt.rcParams['svg.fonttype'] = 'none'
upset(a, sort_by='cardinality', sort_categories_by='cardinality')
plt.savefig(snakemake.output.upset, bbox_inches='tight', dpi=600)
