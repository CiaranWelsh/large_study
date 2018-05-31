# -*- coding: utf8 -*-
from parse import *
import glob, os, pandas
import matplotlib.pyplot as plt
import seaborn
import plotly.plotly as py
import plotly.graph_objs as go
import plotly
import logging
from cycler import cycler
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import matplotlib.ticker as plticker

logging.basicConfig()
LOG = logging.getLogger()
## filename ordered by WG.plate
# design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'

#filenames ordered by batch
# design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design_filename_ordered_by_batch.csv'

# exp = Experiment(design_file)

plotly.tools.set_credentials_file(
    username='c.welsh2',
    api_key='EBoqRg1WvUu2kqsRw1wc'
)




def plot1(treatment='Control', cell_line='B', repeat=list(range(1, 7)), gene='SMAD7', normed=True):
    """

    :param treatment:
    :param cell_line:
    :param repeat:
    :param gene:
    :param normed:
    :return:
    """
    if normed:
        data = exp.treatment_data.loc[treatment, cell_line, repeat, gene]['dct']
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
        treatment_data = treatment_data['dct']
    else:
        treatment_data = treatment_data['Ct']

    treatment_data = treatment_data.swaplevel(1, 2)
    treatment_data = treatment_data.swaplevel(0, 1, axis=1)
    treatment_data = treatment_data.sort_index(level=[0, 1], axis=1)

    # print treatment_data.head()


    for cell_line in exp.cell_lines:
        cell_line_dir = os.path.join(graphs_dir, cell_line)
        for gene in exp.genes:
            print(cell_line, gene)
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

    print(treatment_data)

    for cell in exp.cell_lines:
        for treatment in ['Control', 'TGFb']:
            for gene in exp.genes:
                pass


def plot_individuals(exp, normed=True):
    """

    :param exp:
    :param normed:
    :return:
    """
    if normed:
        data = exp.treatment_data['dct']
    else:
        data = exp.treatment_data['Ct']

    data = data.reorder_levels(['cell_line', 'Assay', 'treatment', 'replicate'])

    graph_dir = os.path.join(exp.root, 'IndividualGraphs')
    for (cell, gene), df in data.groupby(level=[0, 1]):
        cell_dir = os.path.join(graph_dir, cell)
        df = df.loc[cell, gene]
        plt.figure()
        for treat, df2 in df.groupby(level='treatment'):
            treat_dir = os.path.join(cell_dir, treat)
            if not os.path.isdir(treat_dir):
                os.makedirs(treat_dir)
            for rep, df3 in df2.groupby(level='replicate'):
                LOG.debug('plotting {}, {}, {}, {}'.format(cell, gene, treat, rep))
                # print df3.loc[treat, rep]
                plt.rc('axes', prop_cycle=cycler('color', plt.cm.rainbow(numpy.linspace(0, 1, 12))))
                plt.plot(df3.loc[treat, rep].index,
                         df3.loc[treat, rep], label='{}_{}'.format(treat, rep))

                plt.legend(loc=(1, 0.5))
                plt.title('{}_{}'.format(cell, gene))
                if normed:
                    plt.ylabel('2**-(signal - GeoMean(PPIA, B2M))')
                else:
                    plt.ylabel('Ct')

        fname = os.path.join(cell_dir, '{}.png'.format(gene))
        plt.savefig(fname, bbox_inches='tight', dpi=400)
        LOG.debug('plot saved to "{}"'.format(fname))


    def plot_like_preliminary_experiment(data, label, filename):
        """
        Plot the data in the same way as the preliminary data. i.e. normalized to baseline 0
        :param data:
        :param label:
        :param filename:
        :return:
        """
        from matplotlib.backends.backend_pdf import PdfPages
        dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
        design_file = os.path.join(dire, 'new_design.csv')
        #

        E = Experiment(design_file)
        baseline = E.baseline_data.drop('Ct', axis=1)['dct']
        baseline = baseline.drop(96, axis=1)
        baseline.index = baseline.index.swaplevel(0, 1)
        baseline = baseline.loc['Baseline']

        treatment = E.treatment_data.drop('Ct', axis=1)['dct']
        treatment.index = treatment.index.swaplevel('cell_line', 'treatment')

        control = treatment.loc['Control']
        tgfb = treatment.loc['TGFb']

        idx = treatment.index  # .loc['Control']#/baseline
        # control = baseline.loc['Baseline'] / treatment.loc['Control']
        # tgfb = baseline.loc['Baseline'] / treatment.loc['TGFb']
        control_df = {}
        for col in control.columns:
            control_df[col] = control[col] / baseline[0]

        tgfb_dct = {}
        for col in tgfb.columns:
            tgfb_dct[col] = tgfb[col] / baseline[0]

        control_df = pandas.concat(control_df, axis=1)
        tgf_df = pandas.concat(tgfb_dct, axis=1)

        # print control_df.head()
        print(tgf_df.head())

        data = data.reset_index()
        cell_id = []
        for i in data.cell_line:
            if i == 'A' or i == 'B' or i == 'C':
                cell_id.append('Neonatal')

            elif i == 'D' or i == 'E' or i == 'F':
                cell_id.append('Senescent')

            elif i == 'G' or i == 'H' or i == 'I':
                cell_id.append('Adult')

        data['CellType'] = cell_id

        data = data.set_index(['CellType', 'replicate', 'Assay'])
        data = data.drop('cell_line', axis=1)

        data = data.stack().reset_index()
        data = data.rename(columns={'level_3': 'time'})

        sub1_control = control_df.query('cell_line in ["A", "D", "G"]')
        sub1_tgf = tgf_df.query('cell_line in ["A", "D", "G"]')

        sub2_control = control_df.query('cell_line in ["B", "E", "H"]')
        sub2_tgf = tgf_df.query('cell_line in ["B", "E", "H"]')

        sub3_control = control_df.query('cell_line in ["C", "F", "I"]')
        sub3_tgf = tgf_df.query('cell_line in ["C", "F", "I"]')

        dire = r'/home/b3053674/Documents/LargeStudy/GraphsToCompareWithPreliminary'
        dire = r'/home/b3053674/Documents/LargeStudy/Graphs'

        sub1_control_filename = os.path.join(dire, 'ControlSubExp1.pdf')
        sub1_tgf_filename = os.path.join(dire, 'TGFbSubExp1.pdf')
        sub2_control_filename = os.path.join(dire, 'ControlSubExp2.pdf')
        sub2_tgf_filename = os.path.join(dire, 'TGFbSubExp2.pdf')
        sub3_control_filename = os.path.join(dire, 'ControlSubExp3.pdf')
        sub3_tgf_filename = os.path.join(dire, 'TGFbSubExp3.pdf')

        plot_with_time_matched_controls(sub1_control, 'ControlTimeCourseSubexperiment1', sub1_control_filename)
        # plot(sub1_tgf,      'TGFbTimeCourseSubexperiment1:',            sub1_tgf_filename)
        # plot(sub2_control,  'ControlTimeCourseSubexperiment2',          sub2_control_filename)
        # plot(sub2_tgf,      'TGFbTimeCourseSubexperiment2',             sub2_tgf_filename)
        # plot(sub3_control,  'ControlTimeCourseSubexperiment3',          sub3_control_filename)
        # plot(sub3_tgf,      'TGFbTimeCourseSubexperiment3',             sub3_tgf_filename)

        with PdfPages(filename) as pdf:

            for g in E.genes:

                print('plotting ', filename,  g)
                plot_data = data.query('Assay == "{}"'.format(g))
                plot_data = plot_data.rename(columns={0: 'Normed2Ref'})
                if plot_data.empty:
                    print('gene is empty --> ', g)
                    continue

                plt.figure()
                seaborn.barplot(x='time', hue='CellType', y='Normed2Ref', data=plot_data)
                plt.title('{}: {}'.format(label, g))
                plt.legend(loc=(1, 0.3))
                pdf.savefig(bbox_inches='tight', dpi=400)

    def plot_with_time_matched_controls():
        """
        Plot the data in the same way as the preliminary data. i.e. normalized to baseline 0
        :param data:
        :param label:
        :param filename:
        :return:
        """

        from matplotlib.backends.backend_pdf import PdfPages
        dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
        design_file = os.path.join(dire, 'new_design.csv')


        E = Experiment(design_file)
        control = E.treatment_data.query('treatment == "Control"')['dct']
        control.index = control.index.droplevel(1)

        tgfb = E.treatment_data.query('treatment == "TGFb"')['dct']
        tgfb.index = tgfb.index.droplevel(1)


        sub1_control = control.query('cell_line in ["A", "D", "G"]')
        sub1_tgf = tgfb.query('cell_line in ["A", "D", "G"]')

        sub2_control = control.query('cell_line in ["B", "E", "H"]')
        sub2_tgf = tgfb.query('cell_line in ["B", "E", "H"]')

        sub3_control = control.query('cell_line in ["C", "F", "I"]')
        sub3_tgf = tgfb.query('cell_line in ["C", "F", "I"]')

        def plot(ctrl, treat, label, filename, by='time', fold_change=False):
            if fold_change:
                if by == 'time':
                    with PdfPages(filename) as pdf:
                        for g in E.genes:
                            print(label, g)
                            c = ctrl.query('Assay == "{}"'.format(g))
                            t = treat.query('Assay == "{}"'.format(g))
                            if c.empty:
                                continue

                            if t.empty:
                                continue

                            normed = t / c
                            normed = normed.stack().reset_index()

                            if normed.empty:
                                continue

                            cell_id = []
                            for i in normed.cell_line:
                                if i == 'A' or i == 'B' or i == 'C':
                                    cell_id.append('Neonatal')

                                elif i == 'D' or i == 'E' or i == 'F':
                                    cell_id.append('Senescent')

                                elif i == 'G' or i == 'H' or i == 'I':
                                    cell_id.append('Adult')
                            normed['CellType'] = cell_id
                            plt.figure()
                            seaborn.set_context(context='talk', font_scale=1)
                            seaborn.barplot(x='CellType', y=0, hue='time', data=normed)
                            plt.title('{}: {}'.format(label, g))
                            plt.ylabel('Treated / Control (Fold Change)')
                            plt.legend(loc=(1, 0.5))

                            pdf.savefig(bbox_inches='tight', dpi=400)

                elif by =='cell_type':
                    with PdfPages(filename) as pdf:
                        for g in E.genes:
                            print('plotting', label, g)
                            c = ctrl.query('Assay == "{}"'.format(g))
                            t = treat.query('Assay == "{}"'.format(g))

                            if c.empty:
                                continue

                            if t.empty:
                                continue

                            normed = t / c
                            normed = normed.stack().reset_index()

                            if normed.empty:
                                continue

                            cell_id = []


                            for i in normed.cell_line:
                                if i == 'A' or i == 'B' or i == 'C':
                                    cell_id.append('Neonatal')

                                elif i == 'D' or i == 'E' or i == 'F':
                                    cell_id.append('Senescent')

                                elif i == 'G' or i == 'H' or i == 'I':
                                    cell_id.append('Adult')
                            normed['CellType'] = cell_id
                            plt.figure()
                            seaborn.set_context(context='talk', font_scale=1)
                            seaborn.barplot(x='time', y=0, hue='CellType', data=normed)
                            plt.title('{}: {}'.format(label, g))
                            plt.ylabel('Treated / Control (Fold Change)')
                            plt.legend(loc=(1, 0.5))

                            pdf.savefig(bbox_inches='tight', dpi=400)
            if not fold_change:
                if by == 'time':
                    with PdfPages(filename) as pdf:
                        for g in E.genes:
                            print(label, g)
                            c = ctrl.query('Assay == "{}"'.format(g))
                            t = treat.query('Assay == "{}"'.format(g))
                            if c.empty:
                                continue

                            if t.empty:
                                continue

                            t = t.stack().reset_index()

                            if t.empty:
                                continue

                            cell_id = []
                            for i in t.cell_line:
                                if i == 'A' or i == 'B' or i == 'C':
                                    cell_id.append('Neonatal')

                                elif i == 'D' or i == 'E' or i == 'F':
                                    cell_id.append('Senescent')

                                elif i == 'G' or i == 'H' or i == 'I':
                                    cell_id.append('Adult')
                            t['CellType'] = cell_id
                            plt.figure()
                            seaborn.set_context(context='talk', font_scale=1)
                            seaborn.barplot(x='CellType', y=0, hue='time', data=t)
                            plt.title('{}: {}'.format(label, g))
                            plt.ylabel('Treated / Control (Fold Change)')
                            plt.legend(loc=(1, 0.5))

                            pdf.savefig(bbox_inches='tight', dpi=400)

                elif by == 'cell_type':
                    with PdfPages(filename) as pdf:
                        for g in E.genes:
                            print('plotting', label, g)
                            c = ctrl.query('Assay == "{}"'.format(g))
                            t = treat.query('Assay == "{}"'.format(g))

                            if c.empty:
                                continue

                            if t.empty:
                                continue

                            t = t.stack().reset_index()

                            if t.empty:
                                continue

                            cell_id = []

                            for i in t.cell_line:
                                if i == 'A' or i == 'B' or i == 'C':
                                    cell_id.append('Neonatal')

                                elif i == 'D' or i == 'E' or i == 'F':
                                    cell_id.append('Senescent')

                                elif i == 'G' or i == 'H' or i == 'I':
                                    cell_id.append('Adult')
                            t['CellType'] = cell_id
                            plt.figure()
                            seaborn.set_context(context='talk', font_scale=1)
                            seaborn.barplot(x='time', y=0, hue='CellType', data=t)
                            plt.title('{}: {}'.format(label, g))
                            plt.ylabel('Treated / Control (Fold Change)')
                            plt.legend(loc=(1, 0.5))

                            pdf.savefig(bbox_inches='tight', dpi=400)

        dire = '/home/b3053674/Documents/LargeStudy/Graphs/GrpahsNormalizedToTimeMatchedControls'
        sub1_filename_cell_type_fold_change = os.path.join(dire, 'subexperiment1ByCellTypeFoldChange.pdf')
        sub2_filename_cell_type_fold_change = os.path.join(dire, 'subexperiment2ByCellTypeFoldChange.pdf')
        sub3_filename_cell_type_fold_change = os.path.join(dire, 'subexperiment3ByCellTypeFoldChange.pdf')

        sub1_filename_time_fold_change = os.path.join(dire, 'subexperiment1ByTimeFoldChange.pdf')
        sub2_filename_time_fold_change = os.path.join(dire, 'subexperiment2ByTimeFoldChange.pdf')
        sub3_filename_time_fold_change = os.path.join(dire, 'subexperiment3ByTimeFoldChange.pdf')

        sub1_filename_cell_type = os.path.join(dire, 'subexperiment1ByCellType.pdf')
        sub2_filename_cell_type = os.path.join(dire, 'subexperiment2ByCellType.pdf')
        sub3_filename_cell_type = os.path.join(dire, 'subexperiment3ByCellType.pdf')

        sub1_filename_time = os.path.join(dire, 'subexperiment1ByTime.pdf')
        sub2_filename_time = os.path.join(dire, 'subexperiment2ByTime.pdf')
        sub3_filename_time = os.path.join(dire, 'subexperiment3ByTime.pdf')





        plot(sub1_control, sub1_tgf, 'Sub-Experiment 1', sub1_filename_cell_type_fold_change, by='cell_type', fold_change=True)
        plot(sub1_control, sub1_tgf, 'Sub-Experiment 1', sub1_filename_time_fold_change, by='time', fold_change=True)
        plot(sub2_control, sub2_tgf, 'Sub-Experiment 2', sub2_filename_cell_type_fold_change, by='cell_type', fold_change=True)
        plot(sub2_control, sub2_tgf, 'Sub-Experiment 2', sub2_filename_time_fold_change, by='time', fold_change=True)
        plot(sub3_control, sub3_tgf, 'Sub-Experiment 3',sub3_filename_cell_type_fold_change, by='cell_type', fold_change=True)
        plot(sub3_control, sub3_tgf, 'Sub-Experiment 3',sub3_filename_time_fold_change, by='time', fold_change=True)

        plot(sub1_control, sub1_tgf, 'Sub-Experiment 1', sub1_filename_cell_type, by='cell_type', fold_change=False)
        plot(sub1_control, sub1_tgf, 'Sub-Experiment 1', sub1_filename_time, by='time', fold_change=False)
        plot(sub2_control, sub2_tgf, 'Sub-Experiment 2', sub2_filename_cell_type, by='cell_type', fold_change=False)
        plot(sub2_control, sub2_tgf, 'Sub-Experiment 2', sub2_filename_time, by='time', fold_change=False)
        plot(sub3_control, sub3_tgf, 'Sub-Experiment 3',sub3_filename_cell_type, by='cell_type', fold_change=False)
        plot(sub3_control, sub3_tgf, 'Sub-Experiment 3',sub3_filename_time, by='time', fold_change=False)




def plot_baseline():
    pandas.set_option('display.max_rows', 10000)
    seaborn.set(context='talk')
    seaborn.set_style('white')
    seaborn.set_palette(seaborn.color_palette('hls', 2))
    dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
    design_file = os.path.join(dire, 'new_design.csv')

    baseline_dir = r'/home/b3053674/Documents/LargeStudy/Limma/Contrasts/BaselineAnalysis'

    E = Experiment(design_file)
    baseline = E.treatment_data.query('treatment == "Baseline"')['dct']
    baseline = E.baseline_data['dct']
    baseline.index = baseline.index.droplevel(1)
    baseline = pandas.DataFrame(baseline.stack())

    baseline.rename(columns={0: 'DeltaCt'}, inplace=True)
    baseline = baseline.reset_index()

    fname=os.path.join(baseline_dir, 'BaselineGraphs.pdf')
    with PdfPages(fname, 'w') as pdf:
        for g in E.genes[:2]:
            print(g)
            plot_data = baseline.query('Assay == "{}"'.format(g))
            fig = plt.figure()
            ax = seaborn.barplot(x='cell_line', y='DeltaCt', data=plot_data, hue='time',)
            plt.title('Baseline: {}'.format(g))
            plt.ylabel('DeltaCt (mean n=6)')
            plt.xlabel('')
            trans = ax.get_xaxis_transform()
            ax.annotate('Neonatal', xy=(0.4, -0.1), xycoords=trans)

            ax.annotate('', xy=(0.05, -0.06), xycoords='axes fraction', xytext=(0.31, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            ax.annotate('Senescent', xy=(3.5, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.35, -0.06), xycoords='axes fraction', xytext=(0.65, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            ax.annotate('Adult', xy=(6.5, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.68, -0.06), xycoords='axes fraction', xytext=(0.98, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            seaborn.despine(fig=fig, top=False, right=False)
            plt.legend().set_title('Time(h)')
            # pdf.savefig(bbox_inches='tight', dpi=400)
    # print 'file is at "{}"'.format(fname)

def plot_baseline_with_norm_to_neonatal():
    """
    Like above but with Sen and Adult data sets divided by
    neonatal
    :return:
    """
    pandas.set_option('display.max_rows', 10000)
    seaborn.set(context='talk')
    seaborn.set_style('white')
    seaborn.set_palette(seaborn.color_palette('hls', 2))
    dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
    design_file = os.path.join(dire, 'new_design.csv')

    baseline_dir = r'/home/b3053674/Documents/LargeStudy/RecentDataGraphs/Baseline'

    E = Experiment(design_file)
    baseline = E.treatment_data.query('treatment == "Baseline"')['dct']
    baseline = E.baseline_data['dct']
    baseline.index = baseline.index.droplevel(1)
    baseline = pandas.DataFrame(baseline.stack())

    baseline.rename(columns={0: 'DeltaCt'}, inplace=True)

    d = baseline.loc['D'] / baseline.loc['A']
    e = baseline.loc['E'] / baseline.loc['B']
    f = baseline.loc['F'] / baseline.loc['C']
    g = baseline.loc['G'] / baseline.loc['D']
    h = baseline.loc['H'] / baseline.loc['E']
    i = baseline.loc['I'] / baseline.loc['F']

    d['cell_line'] = 'Sen1'
    e['cell_line'] = 'Sen2'
    f['cell_line'] = 'Sen3'
    g['cell_line'] = 'A1'
    h['cell_line'] = 'A2'
    i['cell_line'] = 'A3'
    baseline = pandas.concat([d, e, f, g, h, i])
    # print baseline
    # print d.head()
    baseline = baseline.reset_index()
    fname=os.path.join(baseline_dir, 'BaselineGraphs.pdf')
    with PdfPages(fname, 'w') as pdf:
        for g in E.genes:
            print('plotting', g)
            plot_data = baseline.query('Assay == "{}"'.format(g))
            fig = plt.figure()
            ax = seaborn.barplot(x='cell_line', y='DeltaCt', data=plot_data, hue='time',)
            plt.title('Mean Baseline Levels: {} (n=6)'.format(g))
            plt.ylabel(r'$2^{-\Delta \Delta C_T}$ (err=std)')
            plt.xlabel('')
            trans = ax.get_xaxis_transform()
            # ax.annotate('Neonatal', xy=(0.4, -0.1), xycoords=trans)
            #
            # ax.annotate('', xy=(0.05, -0.06), xycoords='axes fraction', xytext=(0.31, -0.06),
            #             arrowprops=dict(arrowstyle='-',
            #                             color='black',
            #                             linewidth=3)
            #             )

            ax.annotate('Senescent', xy=(1.0, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.05, -0.06), xycoords='axes fraction', xytext=(0.45, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            ax.annotate('Adult', xy=(3.5, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.55, -0.06), xycoords='axes fraction', xytext=(0.95, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            seaborn.despine(fig=fig, top=False, right=False)
            plt.legend().set_title('Time(h)')
            pdf.savefig(bbox_inches='tight', dpi=400)
    print('file is at "{}"'.format(fname))

def plot_baseline_with_norm96_to_0(exp):
    """
    Here instead the 96h time points were calibrated to the 0h time points instead
    :return:
    """
    pandas.set_option('display.max_rows', 10000)
    seaborn.set(context='talk')
    seaborn.set_style('white')
    # seaborn.set_palette(seaborn.color_palette('hls', 2))
    # dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
    # design_file = os.path.join(dire, 'new_design.csv')
    #
    # baseline_dir = r'/home/b3053674/Documents/LargeStudy/RecentDataGraphs/Baseline'
    #
    dct = exp.baseline_data['dct']
    dct.index = dct.index.droplevel('treatment')
    ddct = dct[96] / dct[0]
    baseline_dir = r'/home/b3053674/Documents/LargeStudy/RecentDataGraphs/Baseline'
    ddct = pandas.DataFrame(ddct, columns=['ddct96vs0'])
    ddct = ddct.reset_index()

    print(ddct.head())

    fname=os.path.join(baseline_dir, 'BaselineGraphs96vs9.pdf')
    with PdfPages(fname, 'w') as pdf:
        for g in exp.genes:
            print('plotting', g)
            plot_data = ddct.query('Assay == "{}"'.format(g))
            fig = plt.figure()
            ax = seaborn.barplot(x='cell_line', y='ddct96vs0', data=plot_data)
            plt.title('{} (n=6)'.format(g))
            plt.ylabel(r'$2^{-\Delta \Delta C_T}$ (err=std)')
            plt.xlabel('')
            trans = ax.get_xaxis_transform()
            ax.annotate('Neonatal', xy=(0.4, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.05, -0.06), xycoords='axes fraction', xytext=(0.31, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            ax.annotate('Senescent', xy=(3.5, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.35, -0.06), xycoords='axes fraction', xytext=(0.65, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            ax.annotate('Adult', xy=(6.5, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.68, -0.06), xycoords='axes fraction', xytext=(0.98, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )
            # plt.rc("axes.spines", top=False, right=False)
            seaborn.despine(ax=ax, top=False, right=False)

            pdf.savefig(bbox_inches='tight', dpi=400)
    print('file is at "{}"'.format(fname))


def plot(data, fname):
    """
    must be long form of dataframe without multiindex.
    :param data:
    :param fname:
    :return:
    """
    with PdfPages(fname) as pdf:
        for g in E.genes:
            print(g)
            plot_data = data.query('Assay == "{}"'.format(g))
            print(plot_data)
            fig = plt.figure()
            ax = seaborn.barplot(x='cell_line', y='DeltaCt', data=plot_data, hue='time',)
            plt.title('Baseline: {}'.format(g))
            plt.ylabel('DeltaCt (mean n=6)')
            plt.xlabel('')
            trans = ax.get_xaxis_transform()
            ax.annotate('Neonatal', xy=(0.4, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.05, -0.06), xycoords='axes fraction', xytext=(0.31, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            ax.annotate('Senescent', xy=(3.5, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.35, -0.06), xycoords='axes fraction', xytext=(0.65, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            ax.annotate('Adult', xy=(6.5, -0.1), xycoords=trans)
            ax.annotate('', xy=(0.68, -0.06), xycoords='axes fraction', xytext=(0.98, -0.06),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            seaborn.despine(fig=fig, top=False, right=False)
            plt.legend().set_title('Time(h)')


def plot_ddct(fname=None, baseline_time=0):
    """
    ddct is both control and treated samples
    normalized to the baseline 0 or 96 samples

    :param data:
    :return:
    """
    # pandas.set_option('display.max_rows', 10000)
    from matplotlib.backends.backend_pdf import PdfPages
    seaborn.set(context='talk')
    seaborn.set_style('white')
    seaborn.set_palette(seaborn.color_palette('hls', 2))
    dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
    design_file = os.path.join(dire, 'new_design.csv')

    # baseline_dir = r'/home/b3053674/Documents/LargeStudy/Limma/Contrasts/BaselineAnalysis'
    dire = r'/home/b3053674/Documents/LargeStudy/RecentDataGraphs'
    if fname == None:
        fname = os.path.join(dire, 'ddct_data_{}.pdf'.format(baseline_time))

    E = Experiment(design_file, ddct_type=baseline_time)
    data= deepcopy(E.ddct)
    print(data)
    new_col = []
    for i in range(data.shape[0]):
        new_col.append("{} {}".format(
            str(data.iloc[i]['cell_line']),
            str(data.iloc[i]['treatment'])[0]
        )
        )
    data['x_col'] = new_col
    print(data)
    data = data.sort_values(by='x_col')
    with PdfPages(fname, 'w') as pdf:
        for g in E.genes[:4]:

            try:

                fig, ax, = plt.figure()
                seaborn.despine(fig=fig, top=False, right=False)
                print('plotting', g)
                plot_data = data.query('Assay == "{}"'.format(g))
                new_col = []

                seaborn.barplot(data=plot_data, x='x_col', y=0,
                                hue='time', units='replicate', errwidth=1, ax=ax)
                # plt.xticks(rotation=90)
                plt.title(g)
                plt.xlabel('')
                plt.ylabel('$2^{-\Delta \Delta C_T}$ (err=std)')
                plt.legend(loc=(1, 0.5))
                trans = ax.get_xaxis_transform()
                ax.annotate('Neonatal', xy=(1.4, -0.1), xycoords=trans)
                ax.annotate('', xy=(0.01, -0.06), xycoords='axes fraction', xytext=(0.31, -0.06),
                            arrowprops=dict(arrowstyle='-',
                                            color='black',
                                            linewidth=3)
                            )

                ax.annotate('Senescent', xy=(7.5, -0.1), xycoords=trans)
                ax.annotate('', xy=(0.34, -0.06), xycoords='axes fraction', xytext=(0.65, -0.06),
                            arrowprops=dict(arrowstyle='-',
                                            color='black',
                                            linewidth=3)
                            )

                ax.annotate('Adult', xy=(13.5, -0.1), xycoords=trans)
                ax.annotate('', xy=(0.68, -0.06), xycoords='axes fraction', xytext=(0.98, -0.06),
                            arrowprops=dict(arrowstyle='-',
                                            color='black',
                                            linewidth=3)
                            )

                pdf.savefig(dpi=400, bbox_index='tight')
            except:
                continue
    print('saved to: "{}"'.format(fname))

    ##this is the code I used for running this function:
    # d = r'/home/b3053674/Documents/LargeStudy'
    # fname0 = os.path.join(d, 'ddct_by_baseline0.pdf')
    # fname96 = os.path.join(d, 'ddct_by_baseline96.pdf')
    # fname_av = os.path.join(d, 'ddct_by_baselineAverage.pdf')
    # os.chdir(d)
    # # plot_ddct(baseline_time=0, fname=fname0)
    # plot_ddct(baseline_time=96, fname=fname96)


    # pl(E.control_ddct.reset_index(), 'control_data.pdf')

def plot_d3ct(exp, fname=None, hue='cell_line'):
    """
    calculate d3ct - treated divided by control time course
    :return:
    """
    seaborn.set(context='talk')
    seaborn.set_style('white')
    if fname == None:
        dire = os.path.join(exp.pdf_graph_directory, 'd3ct_data')
        if not os.path.isdir(dire):
            os.makedirs(dire)

        if hue is 'cell_line':
            fname = os.path.join(dire, 'd3ct_data_by_time.pdf')

        elif hue is 'time':
            fname = os.path.join(dire, 'd3ct_data_by_cell_line.pdf')

    d3ct = exp.d3ct
    # data = d3ct.pivot_table(d3ct, index=['cell_line', 'replicate', 'Assay'], columns='time')
    with PdfPages(fname, 'w') as pdf:
        for g in exp.genes[:4]:
            print('plotting {}'.format(g))
            df = d3ct.query('Assay == "{}"'.format(g))
            fig, ax = plt.subplots()
            df = df.reset_index()
            if hue is 'cell_line':
                seaborn.barplot(data=df, x='time', y=0, hue='cell_line', errwidth=1, ax=ax, palette='cubehelix')
                plt.xlabel('Time(h)')
                plt.legend(title='Cell Line')


            elif hue is 'time':
                seaborn.barplot(ax=ax, data=df, x='cell_line', y=0, hue='time', errwidth=1, palette='cubehelix')
                plt.xlabel('')
                plt.legend(title='Time(h)')
                trans = ax.get_xaxis_transform()
                ax.annotate('Neonatal', xy=(0.4, -0.1), xycoords=trans)
                ax.annotate('', xy=(0.05, -0.06), xycoords='axes fraction', xytext=(0.31, -0.06),
                            arrowprops=dict(arrowstyle='-',
                                            color='black',
                                            linewidth=3)
                            )

                ax.annotate('Senescent', xy=(3.5, -0.1), xycoords=trans)
                ax.annotate('', xy=(0.35, -0.06), xycoords='axes fraction', xytext=(0.65, -0.06),
                            arrowprops=dict(arrowstyle='-',
                                            color='black',
                                            linewidth=3)
                            )

                ax.annotate('Adult', xy=(6.5, -0.1), xycoords=trans)
                ax.annotate('', xy=(0.68, -0.06), xycoords='axes fraction', xytext=(0.98, -0.06),
                            arrowprops=dict(arrowstyle='-',
                                            color='black',
                                            linewidth=3)
                            )

            seaborn.despine(fig=fig, ax=ax, top=False, right=False)
            plt.ylabel('TGFb / Control (FC)')
            plt.title(g)
            plt.legend(loc=(1, 0.5))

            pdf.savefig(dpi=400, bbox_inches='tight')
    print('file saved to {}'.format(fname))

def plot_d3ct_bar2(exp, dire, legend=True):
    """
    calculate d3ct - treated divided by control time course
    :return:
    """
    seaborn.set(context='talk', font_scale=2)
    seaborn.set_style('white')

    d3ct = exp.d3ct
    for g in exp.genes:
        print('plotting {}'.format(g))
        df = d3ct.query('Assay == "{}"'.format(g))
        df = df.reset_index()
        fig, ax = plt.subplots()
        pal = seaborn.light_palette('purple', 12, reverse=False)
        seaborn.barplot(data=df, x='cell_line', y=0, hue='time',
                        errwidth=1, palette=pal, linewidth=2, edgecolor='black')
        seaborn.despine(top=True, right=True)
        plt.xlabel('')
        ax.legend_.remove()
        if legend:
            plt.legend(title='Time(h)', loc=(1, -0.1))
        plt.ylabel(r'$\frac{2^{-\Delta\Delta C_{T,TGF\beta}}}{2^{-\Delta\Delta C_{T,Control}}}$')
        plt.title(g)

        trans = ax.get_xaxis_transform()
        ax.annotate('Neonatal', xy=(0.26, -0.17), xycoords=trans)
        ax.annotate('', xy=(0.02, -0.1), xycoords='axes fraction', xytext=(0.32, -0.1),
                    arrowprops=dict(arrowstyle='-',
                                    color='black',
                                    linewidth=3)
                    )

        ax.annotate('Senescent', xy=(2.9, -0.17), xycoords=trans)
        ax.annotate('', xy=(0.34, -0.1), xycoords='axes fraction', xytext=(0.64, -0.1),
                    arrowprops=dict(arrowstyle='-',
                                    color='black',
                                    linewidth=3)
                    )

        ax.annotate('Adult', xy=(6.5, -0.17), xycoords=trans)
        ax.annotate('', xy=(0.68, -0.1), xycoords='axes fraction', xytext=(0.98, -0.1),
                    arrowprops=dict(arrowstyle='-',
                                    color='black',
                                    linewidth=3)
                    )

        fname = os.path.join(dire, '{}.png'.format(g))
        print(fname)
        fig.savefig(fname, dpi=400, bbox_inches='tight')
    print('file saved to {}'.format(fname))


def plot_d3ct_ts(exp, fname=None, hue='cell_line'):
    """
    calculate d3ct - treated divided by control time course
    :return:
    """
    seaborn.set(context='talk', font_scale=2)
    seaborn.set_style('white')
    if fname == None:
        dire = os.path.join(exp.pdf_graph_directory, 'd3ct_data_ts')
        if not os.path.isdir(dire):
            os.makedirs(dire)

        if hue is 'cell_line':
            fname = os.path.join(dire, 'd3ct_data_by_time.pdf')

        elif hue is 'time':
            fname = os.path.join(dire, 'd3ct_data_by_cell_line.pdf')

    d3ct = exp.d3ct
    # data = d3ct.pivot_table(d3ct, index=['cell_line', 'replicate', 'Assay'], columns='time')
    for g in exp.genes:
        print(g)
        try:
            df = d3ct.query('Assay == "{}"'.format(g))
            fig = plt.figure()
            df = df.reset_index()
            df['log2time'] = numpy.log2(df.time)
            ax = seaborn.tsplot(data=df, unit='replicate', time='log2time', condition='cell_line', value=0, color='bright')
                                # err_palette=seaborn.color_palette("hls", 9))

            plt.title(g)
            plt.ylabel(r'$\frac{2^{-\Delta\Delta C_{T,TGF\beta}}}{2^{-\Delta\Delta C_{T,Control}}}$')
            plt.xlabel('Time (log$_2$[h])')
            seaborn.despine(fig=fig, top=True, right=True)
            plt.legend(loc=(1, 0.1))
            fname = os.path.join(dire, '{}.jpg'.format(g))
            print(fname)
            plt.savefig(fname, dpi=400, bbox_inches='tight')
        except Exception as e:
            print(e, e.message)
            continue

def plot_d3ct_ts2(exp):
    """
    Same as plot_d3ct_ts but I am trying to split the time courses up
    so that we have 9 next to each other
    :return:
    """
    seaborn.set(context='paper', font_scale=1)
    seaborn.set_style('white')
    dire = os.path.join(exp.pdf_graph_directory, 'd3ct_data_ts_subplots')
    if not os.path.isdir(dire):
        os.makedirs(dire)

    d3ct = exp.d3ct
    # data = d3ct.pivot_table(d3ct, index=['cell_line', 'replicate', 'Assay'], columns='time')
    for g in exp.genes[:2]:
        print(g)
        # try:
        df = d3ct.query('Assay == "{}"'.format(g))
        fig, ax = plt.subplots(3, 3, sharex='row', sharey='col')
        ax = ax.ravel()
        df = df.reset_index()
        # df['log2time'] = numpy.log2(df.time)
        # df['log10time'] = numpy.log10(df.time)
        cell_lines = df.cell_line.unique()
        bright = [(0.0, 0.24705882352941178, 1.0),
                  (0.0, 0.24705882352941178, 1.0),
                  (0.0, 0.24705882352941178, 1.0),
                  (0.011764705882352941, 0.9294117647058824, 0.22745098039215686),
                  (0.011764705882352941, 0.9294117647058824, 0.22745098039215686),
                  (0.011764705882352941, 0.9294117647058824, 0.22745098039215686),
                  (0.9098039215686274, 0.0, 0.043137254901960784),
                  (0.9098039215686274, 0.0, 0.043137254901960784),
                  (0.9098039215686274, 0.0, 0.043137254901960784)]

        for i in range(len(cell_lines)):
            df2 = df[df.cell_line == cell_lines[i]]
            # plt.subplot(3, 3, i+1)
            # df2.drop('cell_line', inplace=True, axis=1)
            seaborn.tsplot(data=df2, unit='replicate', time='time', value=0, condition='cell_line', color=bright[i],
                           legend=False, ax=ax[i], interpolate=True)
            ax[i].set_title(cell_lines[i], y=0.8, fontsize=16)
            ax[i].set_ylabel('')
            ax[i].set_xlabel('')


            plt.setp(ax[i].xaxis.get_majorticklabels(), rotation=90)
            # loc = plticker.MultipleLocator(base=0.2)
            # ax[i].xaxis.set_major_locator(loc)
            # ax[i].tick_params(axis='both', which='major', labelsize=11)
            # labels = [ for j in ax[i].get_xticklabels()]
            # ax[i].set_xticklabels([0.5, 1, 2, 3, 4, 8, 12, 24, 48, 72, 96])

            # for xpoint, ypoint in zip(x, y):
            #     ax.annotate('{:.2f}'.format(ypoint), (xpoint, ypoint), ha='center',
            #                 va='center', bbox=dict(fc='white', ec='none'))

            if i not in [6, 7, 8]:
                ax[i].get_xaxis().set_visible(False)

            if i not in [0, 3, 6]:
                ax[i].get_yaxis().set_visible(False)

            seaborn.despine(ax=ax[i], top=True, right=True)

        plt.suptitle(g, y=0.95, x=0.4, ha='center', fontsize=20)
        fig.text(0.4, -0.01, 'Time (log$_{10}$[h])', ha='center', fontsize=18)
        fig.text(-0.24, 0.5, r'$\frac{2^{-\Delta\Delta C_{T,TGF\beta}}}{2^{-\Delta\Delta C_{T,Control}}}$', va='center',
                 rotation='vertical', fontsize=20)
        plt.subplots_adjust(left=-0.1, wspace=0.05, hspace=0.05)

        fig.text(-0.17, 0.25, 'Neonatal', rotation='vertical',fontsize=12)
        fig.text(-0.17, 0.55, 'Senescent', rotation='vertical', fontsize=12)
        fig.text(-0.17, 0.75, 'Adult', rotation='vertical', fontsize=12)
        fname = os.path.join(dire, '{}.png'.format(g))
        print(fname)
        plt.savefig(fname, dpi=400, bbox_inches='tight')

        # except Exception as e:
        #     print e, e.message
        #     continue


def plot_d4ct(exp, fname=None, hue='cell_line'):
    """
    calculate d3ct - treated divided by control time course
    :return:
    """
    seaborn.set(context='talk')
    seaborn.set_style('white')
    d4ct = deepcopy(exp.d4ct)
    if fname == None:
        dire = os.path.join(exp.pdf_graph_directory, 'd4ct_data')
        if not os.path.isdir(dire):
            os.makedirs(dire)
        fname = os.path.join(dire, 'd4ct_data.pdf')
    #
    #     if hue is 'cell_line':
    #         fname = os.path.join(dire, 'd3ct_data_by_time.pdf')
    #
    #     elif hue is 'time':
    #         fname = os.path.join(dire, 'd3ct_data_by_cell_line.pdf')

    d4ct.index = d4ct.index.rename(['cell_line', 'replicate','Assay', 'time'])
    data = d4ct.pivot_table(d4ct, index=['cell_line', 'replicate', 'Assay'], columns='time')
    with PdfPages(fname, 'w') as pdf:
        for g in exp.genes:
            print('plotting {}'.format(g))
            df = d4ct.query('Assay == "{}"'.format(g))
            fig = plt.figure()
            df = df.reset_index()
            print(df)

            ax = seaborn.barplot(data=df, x='cell_line', y=0, hue='time', errwidth=1)
            plt.xlabel('')
            plt.legend(title='Time(h)')

            trans = ax.get_xaxis_transform()

            ax.annotate('Senescent', xy=(0.8, -0.12), xycoords=trans)
            ax.annotate('', xy=(0.05, -0.08), xycoords='axes fraction', xytext=(0.45, -0.08),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            ax.annotate('Adult', xy=(4.0, -0.12), xycoords=trans)
            ax.annotate('', xy=(0.55, -0.08), xycoords='axes fraction', xytext=(0.95, -0.08),
                        arrowprops=dict(arrowstyle='-',
                                        color='black',
                                        linewidth=3)
                        )

            # plt.ylabel('TGFb / Control (FC)')
            plt.ylabel(r'$2^{-\Delta ^4C_T}$')
            plt.title(g)
            plt.legend(loc=(1, 0.5))
            seaborn.despine(fig=fig, top=False, right=False)
            pdf.savefig(dpi=400, bbox_inches='tight')
    print('file saved to {}'.format(fname))

def plot_dct_bar(exp, dire, legend=True):
    """
    calculate d3ct - treated divided by control time course
    :return:
    """
    os.makedirs(dire) if not os.path.isdir(dire) else None

    seaborn.set(context='talk', font_scale=2)
    seaborn.set_style('white')

    # print(exp.treatment_data.head())
    # print(sorted(list(set(exp.treatment_data.index.get_level_values(2)))))
    seaborn.set_context('talk', font_scale=4)
    for g in exp.genes:
        print('plotting {}'.format(g))
        treat = exp.treatment_data.query('Assay == "{}"'.format(g))['dct']
        treat.index = treat.index.droplevel('Assay')
        base = exp.baseline_data.query('Assay == "{}"'.format(g))['dct'][0]
        # print(treat.head())
        baseline0 = pandas.DataFrame(base)
        baseline0.index = baseline0.index.droplevel('treatment')
        baseline0.index = baseline0.index.droplevel('Assay')
        # print(baseline0.head())
        control = treat.query("treatment == 'Control'")
        control.index = control.index.droplevel('treatment')

        tgf = treat.query("treatment == 'TGFb'")
        tgf.index = tgf.index.droplevel('treatment')

        control[0] = baseline0
        tgf[0] = baseline0
        control = control.transpose().sort_index(level=[0, 1]).transpose()
        tgf = tgf.transpose().sort_index(level=[0, 1]).transpose()

        control['treatment'] = 'Control'
        tgf['treatment'] = 'TGFB'
        # print(control.head())
        df = pandas.concat([control, tgf]).reset_index().set_index(['cell_line', 'treatment', 'replicate'])
        df = df.stack()
        df = pandas.DataFrame(df)
        df.columns = [r'y']
        df = df.reset_index()

        new_vec = []
        for i in range(df.shape[0]):
            new_vec.append("{}.{}".format(df.iloc[i]['cell_line'], df.iloc[i]['treatment'][:1]))
            # new_vec.append("{}_{}".format(df.iloc[i]['treatment'][:2], df.iloc[i]['cell_line']))

        df['cell_treat'] = new_vec
        df = df.sort_values(by=['cell_line', 'time']).reset_index(drop=True)
        print(df.head())

        fig, ax = plt.subplots(figsize=(30, 15))
        pal = seaborn.light_palette('purple', 12, reverse=False)
        seaborn.barplot(data=df, x='cell_treat', y='y', hue='time',
                        errwidth=1, palette=pal, linewidth=2, edgecolor='black')
        seaborn.despine(top=True, right=True)
        plt.xticks(rotation=90)
        plt.xlabel('')
        ax.legend_.remove()
        if legend:
            plt.legend(title='Time(h)', loc=(1, -0.1))
        plt.ylabel(r'$\Delta C_{T}$')
        plt.title(g)
        # plt.show()

        trans = ax.get_xaxis_transform()
        ax.annotate('Neonatal', xy=(1.5, -0.23), xycoords=trans)
        ax.annotate('', xy=(0.02, -0.15), xycoords='axes fraction', xytext=(0.32, -0.15),
                    arrowprops=dict(arrowstyle='-',
                                    color='black',
                                    linewidth=5)
                    )

        ax.annotate('Senescent', xy=(6.5, -0.23), xycoords=trans)
        ax.annotate('', xy=(0.34, -0.15), xycoords='axes fraction', xytext=(0.64, -0.15),
                    arrowprops=dict(arrowstyle='-',
                                    color='black',
                                    linewidth=5)
                    )

        ax.annotate('Adult', xy=(13.5, -0.23), xycoords=trans)
        ax.annotate('', xy=(0.68, -0.15), xycoords='axes fraction', xytext=(0.98, -0.15),
                    arrowprops=dict(arrowstyle='-',
                                    color='black',
                                    linewidth=5)
                    )

        fname = os.path.join(dire, '{}.png'.format(g))
        print(fname)
        try:
            fig.savefig(fname, dpi=400, bbox_inches='tight')
        except ValueError:
            continue
    print('file saved to {}'.format(fname))



def plot_dct_bar2(exp, dire, legend=True):
    """
    :return:
    """
    seaborn.set(context='talk', font_scale=2)
    seaborn.set_style('white')

    os.makedirs(dire) if not os.path.isdir(dire) else None
    # print(exp.treatment_data.head())
    # print(sorted(list(set(exp.treatment_data.index.get_level_values(2)))))
    for g in exp.genes:
        print('plotting {}'.format(g))
        treat = exp.treatment_data.query('Assay == "{}"'.format(g))['dct']
        treat.index = treat.index.droplevel('Assay')
        base = exp.baseline_data.query('Assay == "{}"'.format(g))['dct'][0]
        # print(treat.head())
        baseline0 = pandas.DataFrame(base)
        baseline0.index = baseline0.index.droplevel('treatment')
        baseline0.index = baseline0.index.droplevel('Assay')
        # print(baseline0.head())
        control = treat.query("treatment == 'Control'")
        control.index = control.index.droplevel('treatment')

        tgf = treat.query("treatment == 'TGFb'")
        tgf.index = tgf.index.droplevel('treatment')

        control[0] = baseline0
        tgf[0] = baseline0
        control = control.transpose().sort_index(level=[0, 1]).transpose()
        tgf = tgf.transpose().sort_index(level=[0, 1]).transpose()
        df = tgf / control

        df = df.stack()
        df = df.reset_index()

        # new_vec = []
    #     for i in range(df.shape[0]):
    #         new_vec.append("{}_{}".format(df.iloc[i]['cell_line'], df.iloc[i]['treatment'][:2]))
    #         # new_vec.append("{}_{}".format(df.iloc[i]['treatment'][:2], df.iloc[i]['cell_line']))
    #
        # df['cell_treat'] = new_vec
        # df = df.sort_values(by='cell_line').reset_index(drop=True)
        print(df.head())
    #
        fig, ax = plt.subplots()
        pal = seaborn.light_palette('purple', 12, reverse=False)
        seaborn.barplot(data=df, x='cell_line', y=0, hue='time',
                        errwidth=1, palette=pal, linewidth=2, edgecolor='black')
        seaborn.despine(top=True, right=True)
        plt.xticks(rotation=90)
        plt.xlabel('')
        ax.legend_.remove()
        if legend:
            plt.legend(title='Time(h)', loc=(1, -0.1))
        plt.ylabel(r'$\Delta C_{T}$')
        plt.title(g)
        # plt.show()

        trans = ax.get_xaxis_transform()
        ax.annotate('Neonatal', xy=(0.26, -0.17), xycoords=trans)
        ax.annotate('', xy=(0.02, -0.1), xycoords='axes fraction', xytext=(0.32, -0.1),
                    arrowprops=dict(arrowstyle='-',
                                    color='black',
                                    linewidth=3)
                    )

        ax.annotate('Senescent', xy=(2.9, -0.17), xycoords=trans)
        ax.annotate('', xy=(0.34, -0.1), xycoords='axes fraction', xytext=(0.64, -0.1),
                    arrowprops=dict(arrowstyle='-',
                                    color='black',
                                    linewidth=3)
                    )

        ax.annotate('Adult', xy=(6.5, -0.17), xycoords=trans)
        ax.annotate('', xy=(0.68, -0.1), xycoords='axes fraction', xytext=(0.98, -0.1),
                    arrowprops=dict(arrowstyle='-',
                                    color='black',
                                    linewidth=3)
                    )

        fname = os.path.join(dire, '{}.png'.format(g))
        print(fname)
        fig.savefig(fname, dpi=400, bbox_inches='tight')
    print('file saved to {}'.format(fname))


if __name__ == '__main__':
    import pickle
    # plot_individuals(exp)
    d = r'/home/b3053674/Documents/LargeStudy'


    exp_pickle = r'/home/b3053674/Documents/LargeStudy/SavedObjects/Experiment.pickle'

    with open(exp_pickle, 'rb') as f:
        exp = pickle.load(f)

    # directory = '/home/b3053674/Documents/LargeStudy/RecentDataGraphs/dct/TGFbByControl'
    # plot_dct_bar2(exp, directory)
    directory = '/home/b3053674/Documents/LargeStudy/RecentDataGraphs/dct/dct'
    plot_dct_bar(exp, directory)


    # plot_d3ct_ts2(exp)#, r'/home/b3053674/Documents/LargeStudy/RecentDataGraphs/d3ct_data_ts_subplots/24-05-2018')

    # fname = os.path.join('/home/b3053674/Documents/LargeStudy/RecentDataGraphs/d3ct_data/24-05-2018', 'cell_line_hue.pdf')
    # dire = '/home/b3053674/Documents/LargeStudy/RecentDataGraphs/d3ct_data/24-05-2018'
    # dire = '/home/b3053674/Documents/LargeStudy/RecentDataGraphs/d3ct_data/24-05-2018/NoLegend'
    # plot_d3ct_bar2(exp, dire=dire, legend=False)




# df = exp.d3ct.query('cell_line == "B" and Assay in ["COL1A1", "COL1A2"]')
    # df = df.unstack()
    # df.index = df.index.droplevel('cell_line')
    # df = df.swaplevel(1, 0)
    # df = df.sortlevel(['Assay'])
    # df = df.transpose()
    # df.index = df.index.droplevel(0)
    # df = df.transpose()
    # df.columns = df.columns * 60
    # df.to_csv('/home/b3053674/Documents/Models/2018/04_April/TGFbModel/Fit6WithOnlyWGData/PlottingForDaryl/COLNeodata.csv')
    # plot_d3ct_ts2(exp)
    # plot_baseline_with_norm96_to_0(exp)







    # fname0 = os.path.join(d, 'ddct_by_baseline0.pdf')
    # fname96 = os.path.join(d, 'ddct_by_baseline96.pdf')
    # fname_av = os.path.join(d, 'ddct_by_baselineAverage.pdf')
    # os.chdir(d)
    # # plot_ddct(baseline_time=0, fname=fname0)
    # plot_ddct(baseline_time=96, fname=fname96)
    #plot_ddct(baseline_time='average', fname=fname_av)
    # dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
    # design_file = os.path.join(dire, 'new_design.csv')

    #
    # pairs_dct1 = {
    #     'D': 'A',
    #     'G': 'A',
    #     'E': 'B',
    #     'H': 'B',
    #     'F': 'C',
    #     'I': 'C',
    # }
    #
    # pairs_dct2 = {
    #     'D': 'C',
    #     'G': 'C',
    #     'E': 'A',
    #     'H': 'A',
    #     'F': 'B',
    #     'I': 'B',
    # }
    #
    # pairs_dct3 = {
    #     'D': 'B',
    #     'G': 'B',
    #     'E': 'C',
    #     'H': 'C',
    #     'F': 'A',
    #     'I': 'A',
    # }

    # E = Experiment(design_file, pairs_dct=pairs_dct1)
    # E = Experiment(design_file, pairs_dct=pairs_dct2)
    # E = Experiment(design_file, pairs_dct=pairs_dct3)

    # print E.d3ct
    # d = r'/home/b3053674/Documents/LargeStudy/RecentDataGraphs/d4ct_data'
    # fname1 = os.path.join(d, 'd4ct_1.pdf')
    # fname2 = os.path.join(d, 'd4ct_2.pdf')
    # fname3 = os.path.join(d, 'd4ct_3.pdf')
    # plot_d4ct(E, fname=fname3)



    # def get_gene(exp, gene, cell_line='A', plot=False, new_name=None,
    #              filename=None):
    #     """
    #     Get fold change averaged data for gene. This is useful as
    #     it formats data for input into copasi
    #     :param exp:
    #     :param gene:
    #     :param cell_line:
    #     :param plot:
    #     :param new_name:
    #     :param filename:
    #     :return:
    #     """
    #     if gene not in exp.genes:
    #         raise ValueError('{} is not in list of genes "{}"'.format(gene, exp.genes))
    #     df = exp.d3ct.query('cell_line == "{}" and Assay == "{}"'.format(cell_line, gene))
    #     df = pandas.pivot_table(df, values=0, index=['cell_line', 'replicate'], columns='time')
    #     mean = df.mean(axis=0)
    #     if new_name is not None:
    #         mean = pandas.DataFrame(mean, columns=[new_name])
    #     else:
    #         mean = pandas.DataFrame(mean, columns=[gene])
    #
    #     mean = mean.reset_index()
    #     mean['time'] = mean['time']*60
    #     mean.loc[0] = [0, 1]
    #
    #     mean['TGFbOn_indep'] = 1
    #
    #     print mean
    #
    #     if filename is not None:
    #         mean.to_csv(filename, sep='\t', index=False)
    #
    #     if plot:
    #         mean.plot(x='time')
    #         plt.show()
    #
    #     return mean

   #
    # smad7 = smad7.mean(axis=0, level='treatment')
    # smad7.loc['fc'] = smad7.loc['TGFb'] / smad7.loc['Control']
    # smad7.transpose().to_csv('/home/b3053674/Documents/Models/2018/02_Feb/StochasticTGFb/Smad7.csv', sep='\t')
    # # print smad7fc
    # # smad7['sum'] = smad7.sum(axis=1)
    # #
    # # smad7 = smad7.div(smad7['sum'], axis=0)
    # # # smad7 = smad7.drop('sum', axis=1)
    # #
    # import matplotlib.pyplot as plt
    # plt.figure()
    # for i in list(smad7.index):
    #     print i
    #     plt.plot(list(smad7.columns),
    #              smad7.loc[i], label='{}'.format(i))
    #     plt.title('{}'.format(i))
    #
    # plt.show()
    #
    # for i in range(1, 7):
    #     plt.plot(list(smad7fc.columns), smad7fc.loc['A', 'SMAD7', i])
    # plt.show()
























