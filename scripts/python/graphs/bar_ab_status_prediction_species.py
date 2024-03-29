#!/usr/bin/python3
import pandas as pd
import functions as ft

""" Stacked bar chart of all predicted expressed/not_expressed snoRNAs per species"""
species_ordered = ['pan_troglodytes', 'gorilla_gorilla', 'macaca_mulatta',
                    'oryctolagus_cuniculus', 'rattus_norvegicus', 'bos_taurus',
                    'ornithorhynchus_anatinus', 'gallus_gallus', 'xenopus_tropicalis',
                    'danio_rerio']
human_df = pd.read_csv(snakemake.input.human_labels, sep='\t')
mouse_df = pd.read_csv(snakemake.input.mouse_labels, sep='\t')
paths = snakemake.input.dfs
dfs = []
for i, path in enumerate(paths):
    species_name = path.split('/')[-1].split('_predicted_label')[0]
    df = pd.read_csv(path, sep='\t')
    df['species_name'] = species_name
    df = df[['predicted_label', 'species_name']]
    dfs.append(df)

# Concat dfs into 1 df
concat_df = pd.concat(dfs)

# Given a species name list, count the number of criteria in specific col of df
# that was previously filtered using species_name_list in global_col
def count_list_species(initial_df, species_name_list, global_col, criteria, specific_col):
    """
    Create a list of lists using initial_col to split the global list and
    specific_col to create the nested lists.
    """
    df_list = []

    #Sort in acending order the unique values in global_col and create a list of
    # df based on these values
    print(species_name_list)
    for val in species_name_list:
        temp_val = initial_df[initial_df[global_col] == val]
        df_list.append(temp_val)


    l = []
    for i, df in enumerate(df_list):
        temp = []
        for j, temp1 in enumerate(criteria):
            crit = df[df[specific_col] == temp1]
            crit = len(crit)
            temp.append(crit)
        l.append(temp)

    return l


# Generate a bar chart of categorical features with a hue of gene_biotype
counts_per_feature = count_list_species(concat_df, species_ordered, 'species_name',
                    list(snakemake.params.hue_color.keys()),
                    'predicted_label')
# Add human and mouse actual labels for comparison
human_expressed = len(human_df[human_df['abundance_cutoff_2'] == 'expressed'])
human_not_expressed = len(human_df[human_df['abundance_cutoff_2'] == 'not_expressed'])
mouse_expressed = len(mouse_df[mouse_df['abundance_cutoff'] == 'expressed'])
mouse_not_expressed = len(mouse_df[mouse_df['abundance_cutoff'] == 'not_expressed'])
counts_per_feature = [[human_expressed, human_not_expressed]] + [[mouse_expressed, mouse_not_expressed]] + counts_per_feature

# Convert to percent
percent = ft.percent_count(counts_per_feature)


# Get the total number of snoRNAs (for which we found snoRNA type) per species
total_nb_sno = str([sum(l) for l in counts_per_feature])
xtick_labels = ['homo_sapiens', 'mus_musculus'] + species_ordered
xtick_labels = [label.capitalize().replace('_', ' ') for label in xtick_labels]


ft.stacked_bar2(percent, xtick_labels,
                list(snakemake.params.hue_color.keys()), '', '',
                'Proportion of snoRNAs (%)', snakemake.params.hue_color, total_nb_sno,
                snakemake.output.bar)
