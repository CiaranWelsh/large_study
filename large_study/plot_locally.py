# -*- coding: utf8 -*-
from parse import *
import glob, os, pandas
import matplotlib.pyplot as plt
import seaborn
import plotly.plotly as py
import plotly.graph_objs as go
import plotly

## filename ordered by WG.plate
design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'

##filenames ordered by batch
design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design_filename_ordered_by_batch.csv'

exp = Experiment(design_file)

plotly.tools.set_credentials_file(
    username='c.welsh2',
    api_key='EBoqRg1WvUu2kqsRw1wc'
)




def plot1(treatment='Control', cell_line=u'B', repeat=range(1, 7), gene='SMAD7', normed=True):
    """

    :param treatment:
    :param cell_line:
    :param repeat:
    :param gene:
    :param normed:
    :return:
    """
    if normed:
        data = exp.treatment_data.loc[treatment, cell_line, repeat, gene]['Norm2ref']
    else:
        data = exp.treatment_data.loc[treatment, cell_line, repeat, gene]['Ct']
    plot_data = data.groupby(level=[0,1]).agg(numpy.mean)
    # trace = go.Scatter(
    #     x=exp.time,
    #     y=plot_data.iloc[0],
    # )
    #
    # pl = [trace]
    #
    # py.iplot(pl, filename='{}_{}_{}_{}'.format(
    #     treatment, cell_line, 2, gene,
    # ))

def plot2(exp, normed=True):
    graphs_dir = os.path.join(exp.root, 'Graphs')
    treatment_data = exp.treatment_data.groupby(level=['cell_line', 'treatment', 'Assay']).agg(
        [numpy.mean, numpy.std]
    )
    if normed:
        treatment_data = treatment_data['Norm2ref']
    else:
        treatment_data = treatment_data['Ct']

    treatment_data = treatment_data.swaplevel(1, 2)
    treatment_data = treatment_data.swaplevel(0, 1, axis=1)
    treatment_data = treatment_data.sort_index(level=[0, 1], axis=1)

    # print treatment_data.head()


    for cell_line in exp.cell_lines:
        cell_line_dir = os.path.join(graphs_dir, cell_line)
        for gene in exp.genes:
            print cell_line, gene
            if not os.path.isdir(cell_line_dir):
                os.makedirs(cell_line_dir)

            c = treatment_data.loc[cell_line, gene, 'Control']['mean']
            cerr = treatment_data.loc[cell_line, gene, 'Control']['mean']
            t = treatment_data.loc[cell_line, gene, 'TGFb']['mean']
            terr = treatment_data.loc[cell_line, gene, 'TGFb']['mean']


            fig, axs = plt.subplots()
            axs.errorbar(c.index, c, yerr=cerr, label='Control')
            axs.errorbar(t.index, t, yerr=terr, label='TGFb')
            title = '{}'.format(
                gene
            )
            plt.legend(loc='best')
            plt.title(title)

            plt.savefig(os.path.join(cell_line_dir, title+'.png'), bbox_inches='tight', dpi=400)

def plot3(exp, normed=True):
    graphs_dir = os.path.join(exp.root, 'Graphs')

    treatment_data = exp.treatment_data.reorder_levels(['cell_line', 'treatment', 'Assay', 'replicate'])
    treatment_data = treatment_data.sort_index(level=[0, 1], axis=1)
    treatment_data = treatment_data.sort_index(level=[0, 1, 2], axis=0)

    print treatment_data

    for cell in exp.cell_lines:
        for treatment in ['Control', 'TGFb']:
            for gene in exp.genes:
                pass


    # for label, df in exp.treatment_data.groupby(level=['cell_line', 'treatment']):
    #     print df, label
    #     if normed:
    #         treatment_data = df.loc[label]['Norm2ref']
    #     else:
    #         treatment_data = df.loc[label]['Ct']
    #
    # print treatment_data
    # treatment_data = treatment_data.swaplevel(1, 2)
    # treatment_data = treatment_data.swaplevel(0, 1, axis=1)
    # treatment_data = treatment_data.sort_index(level=[0, 1], axis=1)
    #
    # # print treatment_data.head()
    #
    #
    # for cell_line in exp.cell_lines:
    #     cell_line_dir = os.path.join(graphs_dir, cell_line)
    #     for gene in exp.genes[:5]:
    #         print cell_line, gene
    #         if not os.path.isdir(cell_line_dir):
    #             os.makedirs(cell_line_dir)
    #
    #         c = treatment_data.loc[cell_line, gene, 'Control']['mean']
    #         cerr = treatment_data.loc[cell_line, gene, 'Control']['mean']
    #         t = treatment_data.loc[cell_line, gene, 'TGFb']['mean']
    #         terr = treatment_data.loc[cell_line, gene, 'TGFb']['mean']
    #
    #         fig, axs = plt.subplots()
    #         axs.errorbar(c.index, c, yerr=cerr, label='Control')
    #         axs.errorbar(t.index, t, yerr=terr, label='TGFb')
    #         title = '{}'.format(
    #             gene
    #         )
    #         plt.legend(loc='best')
    #         plt.title(title)
    #
    #         plt.savefig(os.path.join(cell_line_dir, title + '.png'), bbox_inches='tight', dpi=400)




# def pl(exp):
#     """
#
#     :param exp:
#     :return:
#     """
#     print exp.time
#     print exp.baseline_time
#     print exp.genes
#     print exp.replicates
#
#     for treat in ['Baseline', 'TGFb', 'Control']:
#         if treat is not 'Baseline':
#             for cell in exp.cell_lines:
#                 print exp.treatment_data.loc[treat, cell].groupby(level='replicate').agg([numpy.mean, numpy.std])
                # for rep in exp.replicates:
                #     for gene in exp.genes:
                #         print exp.treatment_data.loc[treat, cell, rep, gene]
    #     if treat is 'Baseline':
    #         pass
    #
    #     else:
    #         print exp.treatment_data.loc[treat]





if __name__ == '__main__':
    plot2(exp)



















