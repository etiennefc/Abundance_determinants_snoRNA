#!/usr/bin/python3
import pandas as pd
import collections as coll

""" Find the snoRNAs that found in at leat 1 test set (across the 10
    iterations) and those that are never found in the test set."""

all_sno_df = pd.read_csv(snakemake.input.all_sno_df, sep='\t')
test_set_paths = snakemake.input.test_sets
test_sets = []
for path in test_set_paths:
    df = pd.read_csv(path, sep='\t')
    test_sets.append(df)

# Count the number of time a snoRNA is present across the 10 iterations
concat_df = pd.concat(test_sets)
sno_occurence_in_test_sets = {}
for sno_id in list(all_sno_df['gene_id_sno']):
    sno_in_test_sets = len(concat_df[concat_df['gene_id_sno'] == sno_id])
    sno_occurence_in_test_sets[sno_id] = sno_in_test_sets

final_df = pd.DataFrame(sno_occurence_in_test_sets.items(), columns=['gene_id_sno', 'nb_occurences_in_test_sets'])
final_df.to_csv(snakemake.output.sno_presence_test_sets, index=False, sep='\t')
