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

logging.basicConfig()
LOG = logging.getLogger()
## filename ordered by WG.plate
design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'

#filenames ordered by batch
# design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design_filename_ordered_by_batch.csv'

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


def plot_individuals(exp, normed=True):
    """

    :param exp:
    :param normed:
    :return:
    """
    if normed:
        data = exp.treatment_data['Norm2ref']
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
        baseline = E.baseline_data.drop('Ct', axis=1)['Norm2ref']
        baseline = baseline.drop(96, axis=1)
        baseline.index = baseline.index.swaplevel(0, 1)
        baseline = baseline.loc['Baseline']

        treatment = E.treatment_data.drop('Ct', axis=1)['Norm2ref']
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
        print tgf_df.head()

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

                print 'plotting ', filename,  g
                plot_data = data.query('Assay == "{}"'.format(g))
                plot_data = plot_data.rename(columns={0: 'Normed2Ref'})
                if plot_data.empty:
                    print 'gene is empty --> ', g
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
        control = E.treatment_data.query('treatment == "Control"')['Norm2ref']
        control.index = control.index.droplevel(1)

        tgfb = E.treatment_data.query('treatment == "TGFb"')['Norm2ref']
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
                            print label, g
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
                            seaborn.barplot(x=u'CellType', y=0, hue='time', data=normed)
                            plt.title('{}: {}'.format(label, g))
                            plt.ylabel('Treated / Control (Fold Change)')
                            plt.legend(loc=(1, 0.5))

                            pdf.savefig(bbox_inches='tight', dpi=400)

                elif by =='cell_type':
                    with PdfPages(filename) as pdf:
                        for g in E.genes:
                            print 'plotting', label, g
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
                            seaborn.barplot(x=u'time', y=0, hue='CellType', data=normed)
                            plt.title('{}: {}'.format(label, g))
                            plt.ylabel('Treated / Control (Fold Change)')
                            plt.legend(loc=(1, 0.5))

                            pdf.savefig(bbox_inches='tight', dpi=400)
            if not fold_change:
                if by == 'time':
                    with PdfPages(filename) as pdf:
                        for g in E.genes:
                            print label, g
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
                            seaborn.barplot(x=u'CellType', y=0, hue='time', data=t)
                            plt.title('{}: {}'.format(label, g))
                            plt.ylabel('Treated / Control (Fold Change)')
                            plt.legend(loc=(1, 0.5))

                            pdf.savefig(bbox_inches='tight', dpi=400)

                elif by == 'cell_type':
                    with PdfPages(filename) as pdf:
                        for g in E.genes:
                            print 'plotting', label, g
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
                            seaborn.barplot(x=u'time', y=0, hue='CellType', data=t)
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
    from matplotlib.backends.backend_pdf import PdfPages
    seaborn.set(context='talk')
    seaborn.set_style('white')
    seaborn.set_palette(seaborn.color_palette('hls', 2))
    dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
    design_file = os.path.join(dire, 'new_design.csv')

    baseline_dir = r'/home/b3053674/Documents/LargeStudy/Limma/Contrasts/BaselineAnalysis'

    E = Experiment(design_file)
    baseline = E.treatment_data.query('treatment == "Baseline"')['Norm2ref']
    baseline = E.baseline_data['Norm2ref']
    baseline.index = baseline.index.droplevel(1)
    baseline = pandas.DataFrame(baseline.stack())

    baseline.rename(columns={0: 'DeltaCt'}, inplace=True)
    baseline = baseline.reset_index()

    with PdfPages(os.path.join(baseline_dir, 'BaselineGraphs.pdf')) as pdf:
        for g in E.genes:
            print g
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

            pdf.savefig(bbox_inches='tight', dpi=400)



    # sub1_control = control.query('cell_line in ["A", "D", "G"]')
    # sub1_tgf = tgfb.query('cell_line in ["A", "D", "G"]')
    #
    # sub2_control = control.query('cell_line in ["B", "E", "H"]')
    # sub2_tgf = tgfb.query('cell_line in ["B", "E", "H"]')
    #
    # sub3_control = control.query('cell_line in ["C", "F", "I"]')
    # sub3_tgf = tgfb.query('cell_line in ["C", "F", "I"]')














if __name__ == '__main__':
    # plot_individuals(exp)

    plot_baseline()



















