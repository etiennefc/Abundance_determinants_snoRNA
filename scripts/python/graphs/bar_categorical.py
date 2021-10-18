#!/usr/bin/python3
import pandas as pd
import functions as ft

df = pd.read_csv(snakemake.input.df, sep='\t')
cd = df[df['sno_type'] == 'C/D']
haca = df[df['sno_type'] == 'H/ACA']

# Generate a bar chart of categorical features with a hue of abundance_cutoff_2
# for all snoRNAs
counts_per_feature = ft.count_list_x(df, 'abundance_cutoff_2',
                    list(snakemake.params.hue_color.keys()),
                    snakemake.wildcards.categorical_features)
percent = ft.percent_count(counts_per_feature)

ft.stacked_bar(percent, sorted(list(df['abundance_cutoff_2'].unique())),
                list(snakemake.params.hue_color.keys()), '', 'Abundance status of snoRNAs',
                'Proportion of snoRNAs (%)', snakemake.params.hue_color,
            snakemake.output.bar_categorical_features)


# Generate a bar chart of categorical features with a hue of abundance_cutoff_2
# for either C/D and H/ACA snoRNAs separately

#For C/D
counts_per_feature_cd = ft.count_list_x(cd, 'abundance_cutoff_2',
                    list(snakemake.params.hue_color.keys()),
                    snakemake.wildcards.categorical_features)
percent_cd = ft.percent_count(counts_per_feature_cd)


ft.stacked_bar(percent_cd, sorted(list(cd['abundance_cutoff_2'].unique())),
                list(snakemake.params.hue_color.keys()), '', 'Abundance status of C/D snoRNAs',
                'Proportion of snoRNAs (%)', snakemake.params.hue_color,
                snakemake.output.bar_categorical_features_cd)


#For H/ACA
counts_per_feature_haca = ft.count_list_x(haca, 'abundance_cutoff_2',
                    list(snakemake.params.hue_color.keys()),
                    snakemake.wildcards.categorical_features)
percent_haca = ft.percent_count(counts_per_feature_haca)

ft.stacked_bar(percent_haca, sorted(list(haca['abundance_cutoff_2'].unique())),
                list(snakemake.params.hue_color.keys()), '', 'Abundance status of H/ACA snoRNAs',
                'Proportion of snoRNAs (%)', snakemake.params.hue_color,
                snakemake.output.bar_categorical_features_haca)
