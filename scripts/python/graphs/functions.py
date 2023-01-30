#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import numpy as np
from sklearn import metrics

"""Functions to create various graphs"""

def heatmap(df, col_color_dict, category, category2, cmap, cbar_label, path, **kwargs):
    """
    Creates a clustered heatmap from a df with a col_color on the side of the graph.
    col_color_dict is a dictionary of the categorical label (keys) and its color
    (values) displayed in the colorbar over the heatmap. category (and category2) are a series
    (of the categorical label) in the SAME ORDER as the data used to create the
    heatmap.
    """

    plt.rcParams['svg.fonttype'] = 'none'

    # Create the two vertical col_colors colorbars over the heatmap
    lut = dict(zip(col_color_dict.keys(), col_color_dict.values()))
    col_colors = category.map(lut)
    col_colors2 = category2.map(lut)

    # Create the heatmap (clustered by rows and columns)
    graph = sns.clustermap(df, cmap=cmap, xticklabels=True,
                            method='weighted', col_colors=pd.concat([col_colors, col_colors2], axis=1), **kwargs)
    plt.xlabel(xlabel=cbar_label, labelpad=-2)
    plt.setp(graph.ax_heatmap.xaxis.get_ticklabels(), rotation=90, fontsize=2)

    # Create legend for colorbar
    legend_list = []
    for i, crit in enumerate(col_color_dict.keys()):
        legend_element = mpatches.Patch(color=list(col_color_dict.values())[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(5, 1),
                fontsize=10)

    plt.savefig(path, bbox_inches='tight', dpi=900)


def heatmap_simple(df, cmap, cbar_label, path, **kwargs):
    """
    Creates a simple clustered heatmap.
    """

    plt.rcParams['svg.fonttype'] = 'none'

    # Create the heatmap (clustered by rows and columns)
    graph = sns.clustermap(df, cmap=cmap, xticklabels=True, yticklabels=True, **kwargs)
    plt.xlabel(xlabel=cbar_label, labelpad=-2)
    #plt.setp(graph.ax_heatmap.xaxis.get_ticklabels(), rotation=90, fontsize=2)


    plt.savefig(path, bbox_inches='tight', dpi=600)


def donut_2(counts, labels, colors, title, legend_labels, legend_colors, path,
            **kwargs):
    """
    Creates a donut plot with two layers. Counts, labels and colors are nested
    lists for outer and inner layers attributes. The counts are represented as
    numerical and percentage values in parentheses.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    ax.axis('equal')
    outer_donut, _, _ = ax.pie(counts[0], radius=1.3, labels=labels[0],
                        colors=colors[0], autopct=make_autopct(counts[0]), pctdistance=0.85, **kwargs)
    plt.setp(outer_donut, width=0.4, edgecolor='white')

    inner_donut, _, _ = ax.pie(counts[1], radius=1.3-0.4, labels=labels[1],
                    colors=colors[1], autopct=make_autopct(counts[1]), pctdistance=0.8, **kwargs)
    plt.setp(inner_donut, width=0.4, edgecolor='white')
    legend_list = []
    for i, crit in enumerate(legend_labels):
        legend_element = mpatches.Patch(color=legend_colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper left', bbox_to_anchor=(-0.1, 1),
                fontsize=10)
    fig.suptitle(title, y=1, fontsize=18)

    plt.savefig(path, dpi=600)


def pie_simple(count_list, colors, title, path, **kwargs):
    """
    Creates a pie chart from a simple list of values.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))

    ax.pie(count_list, colors=list(colors.values()), pctdistance=0.5,
           textprops={'fontsize': 35}, autopct=make_autopct(count_list), **kwargs)
    fig.suptitle(title, y=0.9, fontsize=18)
    plt.legend(labels=colors.keys(), loc='upper right',
                bbox_to_anchor=(1, 1.18), prop={'size': 35})
    plt.savefig(path, dpi=600)

def pie_multiple(count_list, labels, colors, ax_title, title, legend_title, path, **kwargs):
    """
    Creates 8 pie charts from a list of list (count_list) where each global
    element corresponds to a list of local elements (ex: percent per rank across
    8 model intersection).
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(2, 4, figsize=(25, 14))
    plt.subplots_adjust(hspace=0.5)
    ax = axes.flatten()
    for i, element in enumerate(count_list):
        count_per_element = count_list[i][:]
        ax[i].set_title(ax_title[i], fontdict={'fontsize': 25}, x=0.5, y=1.1)
        ax[i].pie(count_per_element, colors=colors, textprops={'fontsize': 19},
                    **kwargs)
        ax[i].axis('equal')
        white_circle = plt.Circle((0, 0), 0.4, color='white') #to create a donut chart
        ax[i].add_artist(white_circle) #to create a donut chart

    fig.suptitle(title, x=0.5, y=0.05, fontsize=22)
    fig.legend(labels=labels, loc='upper right', bbox_to_anchor=(1, 0.5),
                prop={'size': 20}, title=legend_title)
    plt.savefig(path, dpi=600)

def bivariate_density(df, col_x, col_y, hue, path, **kwargs):
    """
    Creates a bivariate density plot (seaborn jointplot) from 2 columns in the df.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    a = sns.jointplot(data=df, x=col_x, y=col_y, hue=hue, kind='scatter', **kwargs)
    plt.savefig(path)


def density(df, xlabel, ylabel, title, path, **kwargs):
    """
    Creates a simple density plot.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.kdeplot(df, shade=True, ax=ax, **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 35})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 35})
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    fig.suptitle(title, fontsize=65, weight='bold', x=0.36, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')


def density_x(df_list, xlabel, ylabel, xscale, title, colors, critère_liste, path, **kwargs):
    """
    Creates a density plot with a number x of dfs to represent on the same ax.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    for i, df in enumerate(df_list):
        sns.kdeplot(df, shade=True, ax=ax, color=colors[i], **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 35})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 35})
    ax.set_xscale(xscale)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)

    legend_list = []
    for i, crit in enumerate(critère_liste):
        legend_element = mpatches.Patch(color=colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0,1.1),
                fontsize=20)

    fig.suptitle(title, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=500)

def density_x_size(df_list, xlabel, ylabel, xscale, title, colors, critère_liste, path, figsize, min, max, **kwargs):
    """
    Creates a density plot with a number x of dfs to represent on the same ax.
    Control the size of the figure and the x-axis limits with the figsize and
    xlim min/max params.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    for i, df in enumerate(df_list):
        sns.kdeplot(df, shade=True, ax=ax, color=colors[i], **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 35})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 35})
    ax.set_xscale(xscale)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    ax.set_xlim(min, max)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)

    legend_list = []
    for i, crit in enumerate(critère_liste):
        legend_element = mpatches.Patch(color=colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper left', bbox_to_anchor=(0,1.1),
                fontsize=30)

    fig.suptitle(title, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=500)

def pairplot(df, hue, hue_color_dict, path, **kwargs):
    """
    Creates a pairplot
    """
    plt.rcParams['svg.fonttype'] = 'none'
    sns.pairplot(df, hue=hue, palette=hue_color_dict)
    plt.savefig(path, bbox_inches='tight', dpi=500)


def scatter(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, optional_annot, path,
            **kwargs):
    """
    Creates a scatter plot (using a x, y and hue column).
    """

    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1,1, figsize=(10,9))
    # ax.set_xscale('symlog')
    #ax.set_yscale('symlog')
    # ax.set_xlim(-0.1, 10000)
    ax.set_ylim(0, 70)
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_title(title, fontsize=10)
    ax.text(3500, 50, optional_annot, fontsize=25)
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col,
                    palette=color_dict, edgecolor='face',
                    s=25, **kwargs)

    plt.legend()
    ax.set_xlabel(xlabel, fontsize=25)
    ax.set_ylabel(ylabel, fontsize=25)
    plt.savefig(path, bbox_inches='tight', dpi=600)


def scatter_cut_axe(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path,
            **kwargs):
    """
    Creates a scatter plot with a y axis cut (using a x, y and hue column).
    """

    plt.rcParams['svg.fonttype'] = 'none'
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(10,9), sharex=True)
    fig.subplots_adjust(hspace=0.05)  # adjust space between axes
    ax1.set_title(title, fontsize=10)

    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col,
                    palette=color_dict, edgecolor='face',
                    s=25, ax=ax1, **kwargs)
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col,
                    palette=color_dict, edgecolor='face',
                    s=25, ax=ax2, **kwargs)
    ax1.set_ylim(64.6, 67)
    ax2.set_ylim(-0.5, 27)
    ax1.tick_params(axis='x', which='both', bottom=False, top=False)
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.legend().remove()
    plt.legend(loc='upper left', bbox_to_anchor=(0, 2), fontsize=20)
    plt.xlabel(xlabel, fontsize=25)
    ax1.set_xlabel(None, fontsize=25)
    ax1.set_ylabel(None, fontsize=25)
    ax2.set_ylabel(None, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)

    # Adding the cut-out slanted lines between the two axes
    suppl_args = dict(marker=[(-1, -0.5), (1, 0.5)], markersize=12,
                linestyle='none', color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0,1], [0,0], transform=ax1.transAxes, **suppl_args)
    ax2.plot([0,1], [1,1], transform=ax2.transAxes, **suppl_args)
    plt.savefig(path, bbox_inches='tight', dpi=600)


def connected_scatter(df, df_hue, hue_col, colors, col_name_for_x_axis,
                        xticklabels, xlabel, ylabel, path, **kwargs):
    """
    Creates a connected scatter plot with multiple lines. The input df must be
    of the following form: a number of lines that correspond to the number of
    gene/model; each line is composed of multiple columns where each column
    corresponds to a categorical value on the x axis (ex: type of dataset (cv,
    train, test)); each line must have a different index (ex: one line per model).
    To create a graph without a hue, you must also create a hue_col but with the
    same value for each line and report only one color.
    """
    # Format the input df so that it is in the good orientation
    # Add also a 'col_name_for_x_axis' column which will correspond to the x axis values
    transposed_df = df.transpose()
    transposed_df[col_name_for_x_axis] = transposed_df.index


    #Colors linked to a potential hue
    print(df_hue[hue_col].unique())
    color_association_dict = dict(zip(list(df_hue[hue_col].unique()), colors))
    print(color_association_dict)
    hue_dict = df_hue.to_dict()
    color_dict = {}
    for key, value in hue_dict.items():
        if key == hue_col:
            for key_2, val_2 in value.items():
                if val_2 in color_association_dict.keys():
                    color_dict[key_2] = color_association_dict[val_2]

    print(hue_dict)
    print(color_dict)


    # Create one line at the time by enumerating across columns
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 30}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.set_xticklabels(xticklabels)

    for i, id in enumerate(transposed_df.iloc[:, :-1].columns):
        sns.lineplot(data=transposed_df, x=col_name_for_x_axis, y=id, marker='o',
                    markeredgecolor=color_dict[id], markerfacecolor=color_dict[id],
                    color=color_dict[id], sort=False, **kwargs)
    ax.set_xlabel(xlabel, fontsize=35)
    ax.set_ylabel(ylabel, fontsize=35)
    ax.set_ylim(0.75, 1)
    legend_list = []
    for i, crit in enumerate(color_association_dict.keys()):
        legend_element = mpatches.Patch(color=color_association_dict[crit], label=crit)
        legend_list.append(legend_element)
    ax.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0.3, 1), fontsize=25)

    plt.savefig(path, bbox_inches='tight', dpi=600)


def connected_scatter_errbars(df, df_hue, hue_col, color_dict, col_name_for_x_axis, stdev,
                        xticklabels, xlabel, ylabel, path, **kwargs):
    """
    Creates a connected scatter plot with multiple lines and error bars. The input
    df must be of the following form: a number of lines that correspond to the
    number of gene/model; each line is composed of multiple columns where each column
    corresponds to a categorical value on the x axis (ex: type of dataset (cv,
    train, test)); each line must have a different index (ex: one line per model).
    To create a graph without a hue, you must also create a hue_col but with the
    same value for each line and report only one color.
    """
    # Format the input df so that it is in the good orientation
    # Add also a 'col_name_for_x_axis' column which will correspond to the x axis values
    transposed_df = df.transpose()
    transposed_df[col_name_for_x_axis] = transposed_df.index

    print(transposed_df)
    print(stdev)
    # Create one line at the time by enumerating across columns
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 30}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.set_xticklabels(xticklabels)

    for i, id in enumerate(transposed_df.iloc[:, :-1].columns):
        sns.lineplot(data=transposed_df, x=col_name_for_x_axis, y=id, marker='o',
                    markeredgecolor=color_dict[id], markerfacecolor=color_dict[id],
                    color=color_dict[id], sort=False, markersize=10, **kwargs)
        ax.errorbar(transposed_df[col_name_for_x_axis], transposed_df[id],
                    [dictio[id] for dictio in stdev], capsize=20, color=color_dict[id])
    ax.set_xlabel(xlabel, fontsize=35)
    ax.set_ylabel(ylabel, fontsize=35)
    ax.set_ylim(0.65, 1.025)
    legend_list = []
    for i, crit in enumerate(color_dict.keys()):
        legend_element = mpatches.Patch(color=color_dict[crit], label=crit)
        legend_list.append(legend_element)
    ax.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0.3, 1), fontsize=25)

    plt.savefig(path, bbox_inches='tight', dpi=600)


def connected_scatter_errbars2(df, df_hue, hue_col, color_dict, col_name_for_x_axis, stdev,
                        xticklabels, xlabel, ylabel, ylim_min, ylim_max, path, **kwargs):
    """
    Creates a connected scatter plot with multiple lines and error bars. The input
    df must be of the following form: a number of lines that correspond to the
    number of gene/model; each line is composed of multiple columns where each column
    corresponds to a categorical value on the x axis (ex: type of dataset (cv,
    train, test)); each line must have a different index (ex: one line per model).
    To create a graph without a hue, you must also create a hue_col but with the
    same value for each line and report only one color.
    """
    # Format the input df so that it is in the good orientation
    # Add also a 'col_name_for_x_axis' column which will correspond to the x axis values
    transposed_df = df.transpose()
    transposed_df[col_name_for_x_axis] = transposed_df.index

    print(transposed_df)
    print(stdev)
    # Create one line at the time by enumerating across columns
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 30}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.set_xticklabels(xticklabels)

    for i, id in enumerate(transposed_df.iloc[:, :-1].columns):
        sns.lineplot(data=transposed_df, x=col_name_for_x_axis, y=id, marker='o',
                    markeredgecolor=color_dict[id], markerfacecolor=color_dict[id],
                    color=color_dict[id], sort=False, markersize=10, **kwargs)
        ax.errorbar(transposed_df[col_name_for_x_axis], transposed_df[id],
                    [dictio[id] for dictio in stdev], capsize=20, color=color_dict[id])
    ax.set_xlabel(xlabel, fontsize=35)
    ax.set_ylabel(ylabel, fontsize=35)
    ax.set_ylim(ylim_min, ylim_max)
    legend_list = []
    for i, crit in enumerate(color_dict.keys()):
        legend_element = mpatches.Patch(color=color_dict[crit], label=crit)
        legend_list.append(legend_element)
    ax.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0.3, 1), fontsize=25)

    plt.savefig(path, bbox_inches='tight', dpi=600)



def violin(df, x_col, y_col, hue_violin, hue_swarm, xlabel, ylabel, title, violin_colors, swarm_colors, path, **kwargs):
    """
    Create a violin plot.
    """
    rc={'ytick.labelsize': 30, 'xtick.labelsize': 25, 'legend.fontsize':25}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    sns.set_context(rc=rc)

    fig, ax = plt.subplots(figsize=(32, 8))
    ax = sns.violinplot(data=df, x=x_col, y=y_col, hue=hue_violin, palette=violin_colors, scale='count', cut=0, **kwargs)
    ax = sns.swarmplot(data=df, x=x_col, y=y_col, hue=hue_swarm, palette=swarm_colors, alpha=0.7, **kwargs)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.title(title)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)

    plt.savefig(path, bbox_inches='tight', dpi=600)


def violin_wo_swarm(df, x_col, y_col, hue_violin, xlabel, ylabel, title, violin_colors, path, **kwargs):
    """
    Create a violin plot.
    """
    rc={'ytick.labelsize': 30, 'xtick.labelsize': 25, 'legend.fontsize':25}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    sns.set_context(rc=rc)

    fig, ax = plt.subplots(figsize=(10, 8))
    ax = sns.violinplot(data=df, x=x_col, y=y_col, hue=hue_violin, palette=violin_colors, scale='count', cut=0, **kwargs)
    plt.title(title)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)

    plt.savefig(path, bbox_inches='tight', dpi=600)


def bar_from_lst_of_lst(list_of_list, bar_positions, bar_colors, width, xticklabels, xlabel, ylabel, hue_list, path, **kwargs):
    """
    Creates a grouped bar chart from a list of lists. Width is generally of 0.3 if the hue is of 2 possibilities) and
    bar_positions is generally equal to [0.15, 0.45] (for 2 bars per tick).
    """
    rc = {'ytick.labelsize': 22}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    for i, type in enumerate(list_of_list):
        x_pos = np.arange(len(type))
        x_pos_i = [x + bar_positions[i] for x in x_pos]
        plt.bar(x_pos_i, type, color=bar_colors[i], width=width, **kwargs)

    #plt.yscale('log')
    plt.xticks([r + width for r in range(len(type))], xticklabels, fontsize=21, rotation=90)
    plt.xlabel(xlabel, fontsize=30)
    plt.ylabel(ylabel, fontsize=30)

    legend_list = []
    for i, crit in enumerate(hue_list):
        legend_element = mpatches.Patch(color=bar_colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper left', bbox_to_anchor=(0, 1), fontsize=25)
    plt.savefig(path, dpi=600, bbox_inches='tight')



def grouped_stacked_bar(lists, x_tick_labels, labels, title, xlabel,
                        ylabel, colors, legend_title, path, **kwargs):
    """
    Create a grouped stacked bar chart. Lists is a list of two lists of lists;
    labels is the labels of the stacked variables.
    """
    plt.rcParams['svg.fonttype'] = 'none'

    # Create dfs from the lists in 'lists'
    df1 = pd.DataFrame(lists[0], index=x_tick_labels, columns=labels)
    df2 = pd.DataFrame(lists[1], index=x_tick_labels, columns=labels)

    # Create the bar plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    df1.plot.bar(ax=ax, position=1, width=.3, color=colors, stacked=True,
                edgecolor='black', **kwargs)
    df2.plot.bar(ax=ax, position=0, width=.3, color=colors, stacked=True,
                edgecolor='black', **kwargs)
    plt.autoscale()
    plt.title(title)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)


    # Add legend
    legend_list = []
    for i, crit in enumerate(labels):
        legend_element = mpatches.Patch(color=colors[crit], label=crit)
        legend_list.append(legend_element)
    legend = ax.legend(handles=legend_list, bbox_to_anchor=(1.1,1.1), fontsize=15)
    legend.set_title(legend_title,prop={'size':18})
    ax.add_artist(legend)

    plt.savefig(path, bbox_inches='tight', dpi=600)


def plot_clustered_stacked2(df_list, colors, xlabel, ylabel, x_tick_labels, path, **kwargs):
    """ Given a list of 2 dataframes (two stacked bars per xtick), create a clustered stacked 
        bar plot. Colors is a list of list of colors for each df. Each df has a column for each variable composing its hue"""

    rc = {'ytick.labelsize': 20, 'xtick.labelsize': 25}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    axe = plt.subplot(111)
    df_list[0].plot(kind='bar', stacked=True, ax=axe, position=-0.05, width=0.35, color=colors[0])
    df_list[1].plot(kind='bar', stacked=True, ax=axe, position=1.05, width=0.35, color=colors[1])
    plt.xlabel(xlabel, fontsize=30)
    plt.ylabel(ylabel, fontsize=30)
    axe.set_xticklabels(x_tick_labels, rotation=0)
    plt.autoscale()
    plt.margins(0.02)

    plt.savefig(path, bbox_inches='tight', dpi=600)


def stacked_bar(lists, x_tick_labels, labels, title, xlabel, ylabel, colors, path, **kwargs):
    """
    Create a stacked bar chart from a list of lists ('lists').
    """
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'

    df = pd.DataFrame(lists, index=x_tick_labels, columns=labels)
    print(df)

    ax = df.plot.bar(stacked=True, figsize=(12,8), color=colors, **kwargs)
    ax.set_xticklabels(x_tick_labels, rotation=0)
    plt.legend(fontsize=35, loc=5, bbox_to_anchor=(0.75, 1.25))
    plt.title(title)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.autoscale()
    plt.savefig(path, bbox_inches='tight', dpi=600)


def stacked_bar2(lists, x_tick_labels, labels, title, xlabel, ylabel, colors, optional_annot, path, **kwargs):
    """
    Create a stacked bar chart from a list of lists ('lists').
    """
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'

    df = pd.DataFrame(lists, index=x_tick_labels, columns=labels)
    print(df)
    ax = df.plot.bar(stacked=True, figsize=(35,8), color=colors, **kwargs)
    ax.set_xticklabels(x_tick_labels, rotation=90)
    # Add optional annotation above bars (ex: number of sno in each bar)
    ax.text(0.5, 105, optional_annot, fontsize=25)
    plt.legend(fontsize=35, loc=5, bbox_to_anchor=(0.75, 1.25))
    plt.title(title)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.autoscale()
    plt.margins(0.02)
    plt.savefig(path, bbox_inches='tight', dpi=600)


def simple_bar(values, x_tick_labels, title, xlabel, ylabel, path, **kwargs):
    """
    Create a simple bar chart from a list of values.
    """
    rc = {'ytick.labelsize': 35, 'xtick.labelsize': 32}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    ax.set_xticks(np.arange(len(values)))
    ax.set_xticklabels(x_tick_labels, rotation=90)
    plt.title(title)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.bar([x for x in range(len(values))], values, **kwargs)
    plt.savefig(path, bbox_inches='tight', dpi=600)


def barh(df, figsize_x, figsize_y, title, xlabel, ylabel, path, **kwargs):
    """
    Creates a stacked horizontal bar chart from a given dataframe.
    """
    rc = {'ytick.labelsize': 15, 'xtick.labelsize': 20}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(figsize_x, figsize_y))
    df.plot.barh(ax=ax, **kwargs)
    plt.legend(fontsize=25, loc=5, bbox_to_anchor=(1, 0.1))
    plt.title(title)
    plt.xlabel(xlabel, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)

    plt.savefig(path, bbox_inches='tight', dpi=600)



def upset_avg_3_cat(sorted_val_vertical_bars, sorted_stdev_vertical_bars,
                    avg_nb, stdev_nb, names, values, vlines_pos, ymins, ymaxs,
                    ylabel_vertical_bars, xlabel_horizontal_bars,
                    ytick_labels_scatter, optional_annot, bar_color, path, **kwargs):
    """
    Create an upset plot with 3 different categories (ex: 3 models) of average
    values (+- stdev) across iterations of these models (not a real upset plot
    per se, since we feed directly the average intersections). The upset plot is
    composed of a vertical bar chart (showing the average size of intersection),
    a scatter plot (under the vertical bar chart, showing what are the
    intersection categories) and a horizontal bar chart showing the average
    number (+- stdev) of items within the 3 different categories.
    **sorted_val_vertical_bars: list of decreasing average values used for the vertical bar chart
    **sorted_stdev_vertical_bars: same as sorted_val_vertical_bars, but standard
                                deviation (across iterations), not the averages
    **avg_nb: list of average values used for the horizontal bar chart (ordered from bottom to top)
    **stdev_nb: same as avg_nb but standard deviation (across iterations), not the averages
    **names: list of all possible x values on the vertical bar chart (ex: TP_TP_TP, TP_FN_TP, ...)
            (the order of the three elements is highly important, as the 1st element
            is for log_reg, the 2nd for svc and the 3rd for rf)
    **values: list of corresponding values (for items in names) on the y axis in
            the scatter plot (ex: y=0 for rf, y=1 for svc and y=2 for log_reg)
    **vlines_pos: list of position of each vertical line on the scatter plot
    **ymins and ymaxs: list of minimal or maximal values for each vertical line
                    on the scatter plot (values are either 0, 2 or 2)
    **ylabel_vertical_bars: ylabel for vertical bar chart
    **xlabel_horizontal_bars: xlabel for horizontal bar chart
    **ytick_labels_scatter: list of y tick labels (from top to bottom) of the scatter
                            plot (ex: ['rf', 'svc', 'log_reg'])
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(2, 2, gridspec_kw={'height_ratios': [2, 1], 'width_ratios': [1, 5]})
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.1, wspace=0.7)
    ax[0, 0].axis('off')

    # Add optional annotation in the first quadrant (ex: "True negatives")
    ax[0, 0].text(0.5, 0.5, optional_annot, fontsize=25)

    # Create vertical bar chart
    ax[0, 1].bar([x for x in range(len(sorted_val_vertical_bars))],
                [val[1] for val in sorted_val_vertical_bars],
                yerr=[val[1] for val in sorted_stdev_vertical_bars], capsize=10,
                color=bar_color)
    ax[0, 1].set_ylabel(ylabel_vertical_bars)

    ## Create scatter with vertical lines linking between vertical dots
    # First, add grey dots for all possible values as background
    ax[1, 1].scatter(list(np.repeat([x for x in range(len(sorted_val_vertical_bars))],
                    len(avg_nb))), [x for x in range(len(avg_nb))] * len(sorted_val_vertical_bars),
                    color='lightgrey', s=75, clip_on=False)
    # Second, add black dots for actual existing values
    ax[1, 1].scatter(names, values, color='black', s=75, clip_on=False)
    # Third, add vertical lines
    ax[1, 1].vlines(vlines_pos, ymins, ymaxs, color='black')

    # Remove/modify ticks and spines for subplots in the upset plot
    ax[0, 1].xaxis.set_tick_params(width=0)
    ax[1, 1].xaxis.set_tick_params(width=0)
    ax[1, 1].yaxis.set_tick_params(width=0)
    ax[1, 0].yaxis.set_tick_params(width=0)
    ax[1, 0].spines['right'].set_color("lightgrey")
    ax[1, 0].spines['left'].set_linewidth(0)
    ax[1, 0].spines['top'].set_linewidth(0)
    ax[0, 1].spines['right'].set_linewidth(0)
    ax[0, 1].spines['bottom'].set_linewidth(0)
    ax[0, 1].spines['top'].set_linewidth(0)
    ax[1, 1].spines['top'].set_linewidth(0)
    ax[1, 1].spines['bottom'].set_linewidth(0)
    ax[1, 1].spines['right'].set_linewidth(0)
    ax[1, 1].spines['left'].set_linewidth(0)

    # Create horizontal bar chart
    ax[1, 0].barh(np.arange(len(avg_nb)), avg_nb, xerr=stdev_nb, capsize=5,
                    color=bar_color, height=0.5)
    ax[1, 0].set_xlabel(xlabel_horizontal_bars)
    ax[1, 0].invert_xaxis()
    ax[1, 0].sharey(ax[1, 1])
    ax[0, 1].sharex(ax[1, 1])
    ax[1, 1].set_xticklabels(['']*len(sorted_val_vertical_bars))
    ytick_labels_scatter.insert(0, "")
    ytick_labels_scatter.append("")
    ax[1, 1].set_yticklabels(ytick_labels_scatter)
    plt.setp(ax[1, 0].get_yticklabels(), visible=False)

    plt.savefig(path, bbox_inches='tight', dpi=600)


def roc_curve(classifier_list, X_test, y_test, xlabel, ylabel, title, path):
    """
    Create a ROC curve from multiple trained classifiers based on their
    performance on the test dataset (X_test, y_test).
    """
    rc = {'ytick.labelsize': 35, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    for i in classifier_list:
        metrics.plot_roc_curve(i, X_test, y_test, ax=ax)
    ax.plot([0.02, 0.98], [0.02, 0.98], transform=ax.transAxes, color='grey', linestyle='--')  # add x=y diagonal
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.title(title, fontsize=40)
    plt.legend(fontsize=30, loc=5, bbox_to_anchor=(1, 0.2))
    plt.margins(x=0.02, y=0.02)
    plt.savefig(path, dpi=600, bbox_inches='tight')


def roc_curve_error_fill(model_name_list, mean_aucs, mean_fpr, mean_tprs, std_tprs,
                        model_colors, output_path, **kwargs):
    """
    Create a roc curve for each average model and add +/- stdev as a cloud
    (error fill) above and below each curve. Each argument is a list where each
    element corresponds to an average model across iterations.
    """
    rc = {'ytick.labelsize': 35, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(15, 15))
    for i, avg_model in enumerate(model_name_list):
        mean_auc = mean_aucs[i]
        mean_auc = str(np.round(mean_auc, 3))
        ax.plot(mean_fpr, mean_tprs[i],
                label=f'{avg_model} (mean AUC = {mean_auc})',
                lw=1.5, color=model_colors[avg_model])
        std_tpr_upper = np.minimum(mean_tprs[i] + std_tprs[i], 1)
        std_tpr_lower = np.maximum(mean_tprs[i] - std_tprs[i], 0)
        ax.fill_between(mean_fpr, std_tpr_upper, std_tpr_lower, color=model_colors[avg_model], alpha=0.2)
    ax.plot([0.02, 0.98], [0.02, 0.98], transform=ax.transAxes, color='grey', linestyle='--')
    plt.margins(x=0.02, y=0.02)
    plt.xlabel("False positive rate", fontsize=40)
    plt.ylabel("True positive rate", fontsize=40)
    ax.legend(loc="lower right", fontsize=30)
    plt.savefig(output_path, dpi=600, bbox_inches='tight')


def count_list_x(initial_df, global_col, criteria, specific_col):
    """
    Create a list of lists using initial_col to split the global list and
    specific_col to create the nested lists.
    """
    df_list = []

    #Sort in acending order the unique values in global_col and create a list of
    # df based on these values
    print(sorted(list(initial_df[global_col].unique())))
    for val in sorted(list(initial_df[global_col].unique())):
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

def count_list_x_unsorted(initial_df, global_label_list, global_col, criteria, specific_col):
    """
    Create a list of lists using initial_col to split the global list and
    specific_col to create the nested lists.
    """
    df_list = []

    #Sort in acending order the unique values in global_col and create a list of
    # df based on these values
    print(global_label_list)
    for val in global_label_list:
        temp_val = initial_df[initial_df[global_col] == val]
        df_list.append(temp_val)


    l = []
    print(criteria)
    for i, df in enumerate(df_list):
        temp = []
        for j, temp1 in enumerate(criteria):
            crit = df[df[specific_col] == temp1]
            crit = len(crit)
            temp.append(crit)
        l.append(temp)

    return l


def percent_count(count_list):
    """
    Create from a list of lists a percentage list of lists.
    """
    percent_list = []

    for i, cat in enumerate(count_list):
        temp = []
        for j, crit in enumerate(cat):
            total = sum(cat)
            if total != 0:
                percent = crit/total * 100
                percent = round(percent, 2)
            else:
                percent = 0
            temp.append(percent)
        percent_list.append(temp)

    return percent_list


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d} \n ({p:.1f}%)'.format(p=pct,v=val)
    return my_autopct

def fisher_contingency(group1, group2, col, crit):
    count1a = len(group1[group1[col] == crit])
    count1b = len(group1[group1[col] != crit])
    count2a = len(group2[group2[col] == crit])
    count2b = len(group2[group2[col] != crit])

    dict = {'group1': [count1a, count1b], 'group2': [count2a, count2b]}
    table = pd.DataFrame(data=dict, index=[crit, '!= '+crit])
    print(table)
    return table

'''
def violin_df(df_list, common_col_str, criteria):
    """
    Create a df suitable for the violin plot function of Seaborn.
    """
    for i, xtick in enumerate(criteria)



    final_df = pd.concat(df_list, axis=0)





    return final_df
'''
