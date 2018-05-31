from .parse import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import Imputer
import matplotlib.pyplot as plt
import seaborn
import logging
from cycler import cycler
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly
from flask import Flask
import pickle



def do_pca_ddct(exp, thresh=3, exclude_query=None):

    ddct = exp.ddct.set_index(['treatment', 'cell_line', 'Assay', 'replicate', 'time'])
    ddct = ddct.rename(columns={0:'ddct'})
    design = exp.design[['Sample', 'Treatment Start Date', 'Filename']]


    ddct = ddct.reset_index() ##treatment cell_line    Assay  replicate  time      ddct
    design = design.reset_index() ##sub_experiment cell_id treatment  replicate  time_point  Sample Treatment Start Date                  Filename

    design = design.rename(columns={'time_point': 'time'})
    ddct = ddct.rename(columns={'cell_line': 'cell_id'})
    idx = ['treatment', 'cell_id', 'replicate', 'time']
    ddct = ddct.set_index(idx)
    design = design.set_index(idx)

    df = design.merge(ddct, left_index=True, right_index=True)
    # print df
    # df.to_csv('/home/b3053674/Documents/LargeStudy/SavedObjects/data.csv')


    idx = idx + ['Filename', 'Assay', 'Treatment Start Date', 'Sample']
    df =df.reset_index().set_index(idx)


    df = df.reset_index()

    cell_line = []
    for i in df['cell_id']:
        if i in ['A', 'B', 'C']:
            cell_line.append('Neonatal')

        elif i in ['D', 'E', 'F']:
            cell_line.append('Senescent')

        elif i in ['G', 'H', 'I']:
            cell_line.append('Adult')

    ## add cell line column
    df['cell_line'] = cell_line

    df = df.pivot_table(columns=['Sample', 'cell_id', 'cell_line', 'replicate',
                                 'treatment', 'time', 'Treatment Start Date',
                                 'sub_experiment', 'Filename'], index='Assay')
    df.dropna(axis=0, how='all', thresh=thresh, inplace=True)

    I = Imputer(axis=1, strategy='median')
    imputed = pandas.DataFrame(I.fit_transform(df), index=df.index, columns=df.columns)
    imputed = imputed.transpose()

    pca = PCA(10)
    pca.fit(imputed)
    explained_var = pandas.DataFrame(pca.explained_variance_ratio_)
    explained_var.to_pickle('/home/b3053674/Documents/LargeStudy/SavedObjects/PCAExplainedVar.pickle')
    pc = pca.transform(imputed)
    pc = pandas.DataFrame(pc)

    # print explained_var

    pc.index = imputed.index
    pc = pc[[0, 1]]

    return pc, explained_var


def do_pca_ct(exp, thresh=3, exclude_query=None):

    treatment = exp.treatment_data.stack()
    baseline = pandas.DataFrame(exp.baseline_data).stack()
    ct = pandas.concat([treatment, baseline])['Ct']
    ct = pandas.DataFrame(ct, columns=['Ct'])

    # ct = ct.reset_index() ##treatment cell_line    Assay  replicate  time      ddct
    design = exp.design[['Sample', 'Treatment Start Date', 'Filename']]

    design = design.reset_index() ##sub_experiment cell_id treatment  replicate  time_point  Sample Treatment Start Date                  Filename

    design = design.rename(columns={'time_point': 'time'})
    ## set index to be the same in both so that we can merge the design with data
    index_cols = ['cell_id', 'treatment', 'replicate', 'time']
    ct = ct.reset_index().rename(columns={'cell_line': 'cell_id'}).set_index(index_cols)
    design = design.reset_index().set_index(index_cols)
    df = design.merge(ct, left_index=True, right_index=True)
    del df['index']

    df = df.reset_index()

    cell_line = []
    for i in df['cell_id']:
        if i in ['A', 'B', 'C']:
            cell_line.append('Neonatal')

        elif i in ['D', 'E', 'F']:
            cell_line.append('Senescent')

        elif i in ['G', 'H', 'I']:
            cell_line.append('Adult')

    ## add cell line column
    df['cell_line'] = cell_line
    df = df.pivot_table(columns=['Sample', 'cell_line', 'cell_id', 'replicate',
                                 'treatment', 'time', 'Treatment Start Date',
                                 'sub_experiment', 'Filename'], index='Assay')

    df.dropna(axis=0, how='all', thresh=thresh, inplace=True)

    I = Imputer(axis=1, strategy='median')
    imputed = pandas.DataFrame(I.fit_transform(df), index=df.index, columns=df.columns)
    imputed = imputed.transpose()
    imputed.index = imputed.index.droplevel(0)

    pca = PCA(10)
    pca.fit(imputed)
    explained_var = pandas.DataFrame(pca.explained_variance_ratio_)

    pc = pca.transform(imputed)
    pc = pandas.DataFrame(pc)

    pc.index = imputed.index
    pc = pc[[0, 1]]
    return pc, explained_var


def plot_pca(pc, explained_var, colour_by='treatment', title=None,
             folder=r'/home/b3053674/Documents/LargeStudy/RecentDataGraphs/QC/PCAPlots',
             query=None):
    """
    Plot PCA
    :param pc:
    :param explained_var:
    :param colour_by:
    :param title:
    :param folder:
    :param query:
    :return:
    """

    if title is None:
        title = colour_by

    if query is not None:
        pc = pc.query(query)

    seaborn.set_context('talk', font_scale=2)
    seaborn.set_style('white')
    fig = plt.figure()


    colours = iter(seaborn.color_palette('hls', len(pc.groupby(colour_by))+2))

    for label, df in pc.groupby(level=colour_by):
        plt.plot(df[0], df[1], 'o', alpha=0.5, color=next(colours), label=label)
        plt.xlabel('PC1 ({:.2f}%)'.format(explained_var[0].iloc[0]*100))
        plt.ylabel('PC2 ({:.2f}%)'.format(explained_var[0].iloc[1]*100))
        plt.title('PCA coloured by {}'.format(colour_by))

    if colour_by is 'Filename':
        plt.legend(loc=(1, 0.1), title=colour_by, fontsize=15)
    else:
        plt.legend(loc=(1, 0.1), title=colour_by, fontsize=25)

    seaborn.despine(fig, top=True, right=True)
    fname = os.path.join(folder, colour_by+'.png')
    print('saved to "{}"'.format(fname))
    plt.savefig(fname, dpi=800, bbox_inches='tight')

def plot_scree(explained_var, folder):
    explained_var = explained_var.reset_index()
    explained_var['index'] = explained_var['index'] + 1
    explained_var[0] = explained_var[0]*100
    print(explained_var)
    seaborn.set_style('white')
    seaborn.set_context('talk', font_scale=2)
    fig = plt.figure()
    seaborn.barplot(data=explained_var, x='index', y=0)
    plt.xlabel('Principle Components (PC)')
    plt.ylabel('Percentage Variance \nExplained by PC')
    seaborn.despine(fig, top=True, right=True)
    fname = os.path.join(folder, 'scree_plot.png')
    plt.savefig(fname, dpi=500, bbox_inches='tight')
    print('saved to "{}"'.format(fname))

if __name__ == '__main__':
    experiment_pickle = '/home/b3053674/Documents/LargeStudy/SavedObjects/Experiment.pickle'

    with open(experiment_pickle, 'rb') as f:
        exp = pickle.load(f)

    factors = ['cell_id', 'cell_line', 'replicate', 'treatment', 'time',
               'Treatment Start Date', 'sub_experiment', 'Filename']
    pc, explained_var = do_pca_ct(exp)
    # pc, explained_var = do_pca_ddct(exp)

    raw_ct_folder = '/home/b3053674/Documents/LargeStudy/RecentDataGraphs/QC/PCAPlots/RawDataQC'
    # ddct_folder = '/home/b3053674/Documents/LargeStudy/RecentDataGraphs/QC/PCAPlots/ddctDataQC'

    # for f in factors:
    #     plot_pca(pc, explained_var, colour_by=f, folder=ddct_folder)
    plot_scree(explained_var, folder=raw_ct_folder)

