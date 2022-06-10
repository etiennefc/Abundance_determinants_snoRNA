#!/usr/bin/python3
import pandas as pd
import functions as ft
from scipy import stats as st

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
total_nb_sno = [sum(l) for l in counts_per_feature]
xtick_labels = ['homo_sapiens', 'mus_musculus'] + species_ordered
sno_nb_dict = dict(zip(xtick_labels, total_nb_sno))

# Create df
df = pd.DataFrame(percent, index=xtick_labels, columns=list(snakemake.params.hue_color.keys()))
df = df.reset_index()
df = df.rename(columns={'index': 'species'})

# Create sno_nb and predicted_vs_actual_ab_status cols
df['sno_nb'] = df['species'].map(sno_nb_dict)
df.loc[(df.species == 'homo_sapiens') | (df.species == 'mus_musculus'), 'predicted_vs_actual_ab_status'] = 'Actual abundance status'
df.predicted_vs_actual_ab_status = df.predicted_vs_actual_ab_status.fillna('Predicted abundance status')
print(df)

# Create scatter plot
color_dictio = {'Actual abundance status': '#000000',
                'Predicted abundance status': '#bdbdbd'}
pearson_r, pval = st.pearsonr(list(df.sno_nb), list(df.expressed))
print(pearson_r, pval)
ft.scatter(df, 'sno_nb', 'expressed', 'predicted_vs_actual_ab_status',
            'Number of snoRNAs per species', 'Proportion of expressed snoRNAs (%)',
            '', color_dictio, f"Pearson's r: {pearson_r}\np-value: {pval}", snakemake.output.scatter)
