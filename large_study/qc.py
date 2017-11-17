from parse import *
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

logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger()


server = Flask(__name__)

app = dash.Dash(server=server)

app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})


design_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'GSS2375_WB_NewDur_Grant')
design_file = os.path.join(design_file, 'new_design.csv')
if not os.path.isfile(design_file):
    os.chdir('..')
    design_file = os.path.abspath(design_file)


if not os.path.isfile(design_file):
    raise ValueError('"{}" not valid file'.format(design_file))

# design_file = os.path
exp = Experiment(design_file)
sample_keys = exp.subexperiments[1].plates[1].samples.keys()
genes = exp.subexperiments[1].plates[1].samples[sample_keys[0]].genes


# app.layout = html.Div([
#
# ])


'''
First do a PCA on the 3 subexperiments, sets of 6 
'''
print exp.subexperiments.plates




















#
#
# seaborn.set_context(context='talk', font_scale=1.5)
#
#
# plt.rc('axes', prop_cycle=cycler('color', ['r', 'g', 'b', 'k', 'y', 'c']))
#
#
#

#
# exp = Experiment(design_file)
#
# qc_dir = os.path.join(exp.root, 'QC')
# if not os.path.isdir(qc_dir):
#     os.makedirs(qc_dir)
#
#
# def pca(exp, genes=None):
#     """
#
#     :param exp:
#     :param genes:
#     :return:
#     """
#     if genes is not None:
#         if not isinstance(genes, list):
#             genes = [genes]
#     data = exp.treatment_data['Norm2ref']#.dropna(axis=1, how='all')#, thresh=5)
#
#     LOG.debug('shape before dropping NaNs --> {}'.format(data.shape))
#     data.dropna(axis=0, how='all', thresh=8, inplace=True)
#     LOG.debug('shape after dropping NaNs --> {}'.format(data.shape))
#     I = Imputer(axis=1)
#     data = pandas.DataFrame(I.fit_transform(data), index=data.index, columns=data.columns)
#     # print data.head()
#     data = data.reorder_levels(['cell_line', 'treatment', 'Assay', 'replicate'])
#
#
#     pc_dict = {}
#     for cell in exp.cell_lines:
#         cell_dir = os.path.join(qc_dir, cell)
#         for treat in ['Control', 'TGFb']:
#             fig=plt.figure()
#             treat_dir = os.path.join(cell_dir, treat)
#             if not os.path.isdir(treat_dir):
#                 os.makedirs(treat_dir)
#             if genes is None:
#                 genes = list(set(data.loc[cell, treat].index.get_level_values(0)))
#             for gene in genes:
#                 gene_dir = os.path.join(treat_dir, gene)
#
#                 # print data.loc[cell, treat, gene]
#                 # print cell, treat, gene
#                 # print data.loc[cell, treat, gene]
#                 try:
#                     pc = PCA(2).fit_transform(data.loc[cell, treat, gene])
#                     pc = pandas.DataFrame(pc, index=data.loc[cell, treat, gene].index)
#                     plt.scatter(pc[0], pc[1], label='{}_{}'.format(gene, treat))
#                     plt.legend(loc=(1, 0.5))
#                 except KeyError:
#                     continue
#                 # try:
#                 #     for i in range(pc.shape[0]):
#                 #         plt.scatter(pc.iloc[i][0], pc.iloc[i][1], marker='o', label=str(pc.iloc[i].name))
#                 #         plt.legend(loc=(0.5, 1))
#                 # except KeyError:
#                 #     print 'key error', cell, treat, gene
#                 # plt.plot(pc[0], pc[1], label='{}_{}_{}'.format(cell, treat, gene))
#                 # for i in range(pc.shape[0]):
#                 #     x = pc.iloc[i][0]
#                 #     y = pc.iloc[i][1]
#                 #     plt.plot(x, y, marker='o', label='{}_{}_{}'.format(cell, treat, gene))
#                 #     plt.legend(loc=(1,1))
#
#
#             fname = os.path.join(treat_dir, reduce(
#                 lambda x, y: '{}_{}'.format(x, y), genes)
#                                  )
#             fname = fname+'.png'
#             plt.savefig(fname, bbox_inches='tight', dpi=400)
#             LOG.debug('saved to "{}"'.format(fname))
#
#
#
# def get_gene_batches(n=6):
#     """
#     PCA plot gets crowded quite quickly. This
#     function splits the 72*1 gene list into
#     6*12 gene list. Then each of the 6 sets of
#     12 will be PCA'ed together
#     :return:
#     """
#     len_genes = len(exp.genes)
#     # if len_genes % n is not 0:
#     #     raise ValueError('{} / {} not exactly divisible. Pick another n number'.format(
#     #         len(exp.genes), n
#     #     ))
#
#     num_per_batch = len_genes / n
#     gene_list = []
#     for i in range(num_per_batch):
#         LOG.debug('i is --> {}'.format(i))
#         gene_list2 = []
#         for j in range(n):
#             LOG.debug('j is --> {}'.format(j))
#             index = i*n+j
#             gene_list2.append(exp.genes[index])
#         gene_list.append(gene_list2)
#     return gene_list
#
#
#


if __name__ == '__main__':
    pass
    # app.run_server(debug=True)
    # # pca(exp, ['SMAD7', 'SMAD3'])
    #
    # gene_list = get_gene_batches()
    # for i in gene_list:
    #     pca(exp, i)
    # # print gene_list
    # print len(gene_list)
    # print len(gene_list[0])









