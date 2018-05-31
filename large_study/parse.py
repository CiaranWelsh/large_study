# -*- coding: utf8 -*-
"""
Analysis of Wafergen data for LargeStudy
========================================

The Experiment
--------------

Fibroblast were treated with 5ng/mL TGFb or control or nothing (baseline).

Fibroblasts were one of 9 groups:

* Neonatal Human Dermal Fibroblasts (HDFn) * 3
* Irradiated Human Dermal Fibroblasts (and therefore Senescent) (IR) * 3
* Adult Human Dermal Fibroblasts (A) * 3

The HDFn and IR cell lines are matched.

Variables:
* Cell line 3*3 as above
* Time Points: 0.5, 1, 2,
* Replicates: 6
* Groups TGFb, control and baseline (latter only has 0 and 96h)


This module contains five classes:

* Experiment
* SubExperiment
* Plate
* Sample
* Query

The Experiment has three sub-experiments, each consisting of a HDFn, IR
and A cell line apiece. Each sub-experiment has 6 plates, each containing
all time points, groups and one replicate of each of the three cell lines. Each
plate has 72 Samples and 72 genes. Each Sample holds the data belonging
to that sample and normalizes to the geometric mean of the reference
genes.

The Query class enables easy searching for the subset of data needed.

The baseline is split from the Control and TGFb treatment groups
because it only has the 0 and 96h time points. Each of the :py:class:`Experiment`,
:py:class:`SubExperiment` and :py:class:`Plate` have two properties:

# treatment_data: Get the control and TGFb data
# baseline_data: Get the baseline data
"""


import pandas, numpy
import seaborn
import matplotlib.pyplot as plt
import os, glob
from copy import deepcopy
from scipy.stats.mstats import gmean
import logging
from cached_property import cached_property
import pickle
from functools import reduce
# FORMAT = "%(asctime)s"
logging.basicConfig(level=logging.DEBUG)

pandas.options.mode.chained_assignment = None



class Experiment(object):
    """
    Class to hold the entire experiment. Every aspect
    of the experiment is obtainable via this class.
    Accepts a design file. The design file had date and
    filename columns manually added before processing.
    """
    def __init__(self, design, max_nan=6, ddct_type=0, from_files=True, pdf_graph_directory=None, pairs_dct=None,
                 outliers=True):
        """

        :param design:
            Path to the design file.

        :param max_nan:
            If max_nan NaN's present in treated or control profile,
            remove the data entry. default = 6.

        :param ddct_type:
        :param from_files:
        :param_pdf_graph_directory:
        :param pairs_dct
        """
        self.pairs_dct = pairs_dct
        self._design = design
        self.from_files = from_files
        self.max_nan = max_nan
        self.root = os.path.dirname(design)
        self.subexperiments = self.create_subexperiments()
        self.ddct_type = ddct_type
        self.pdf_graph_directory = pdf_graph_directory
        self.outliers = outliers

        if self.pdf_graph_directory is None:
            self.pdf_graph_directory = '/home/b3053674/Documents/LargeStudy/RecentDataGraphs'
        if self.ddct_type not in [0, 96, 'average']:
            raise TypeError('ddct_type should be 0, 96 or average')
        ## nested unzipping of nested list comprehension to retrieve all sample names from all plates
        self.all_samples = reduce(
            lambda x, y: x + y,
            [list(self.subexperiments[i].plates[j].samples.keys()) for i in self.subexperiments for j in self.subexperiments[i].plates]
        )

        if len(self.all_samples) != 1296:
            raise Exception('length of all samples should be 1296, not "{}"'.format(len(self.all_samples)))

        self.cell_lines = sorted(reduce(
            lambda x, y: x+y, [self.subexperiments[i].cell_id for i in self.subexperiments]
        ))
        self.genes = sorted(list(set(self.treatment_data.index.get_level_values(3))))
        self.cell_lines = sorted(list(set(self.treatment_data.index.get_level_values(0))))
        self.replicates = sorted(list(set(self.treatment_data.index.get_level_values(2))))

        self.time = sorted(list(set(self.treatment_data.columns.get_level_values(1))))
        self.baseline_time = sorted(list(set(self.baseline_data.columns.get_level_values(1))))
        if self.pairs_dct is None:
            self.pairs_dct = {
                    'D': 'A',
                    'G': 'A',
                    'E': 'B',
                    'H': 'B',
                    'F': 'C',
                    'I': 'C',
                }

        self.ddct = self.calculate_ddct(self.ddct_type)

        ## remove outlers if self.outliers is false
        if not self.outliers:
            print('removeing outliers')
            self.ddct = self.remove_outliers()

        self.d3ct = self.calculate_d3ct()
        self.d4ct = self.calculate_d4ct(self.pairs_dct)

    def __getitem__(self, item):
        return self.subexperiments[item]

    def remove_outliers(self):
        """

        :return:
        """
        def remove(df, pattern):
            """
            Return df with all enteries which df.query(pattern) match with removed
            :param df:
            :param pattern:
            :return:
            """
            return df[~df.index.isin(df.query(pattern).index)]

        ## elastin data is causing cell line G control, 3h to look anomylous on PCA.
        ## After further investigation I have deemed the following conditions unreliable data
        pattern1 = "cell_line in ['G', 'D', 'I'] and " \
                   "treatment == 'Control' and " \
                   "Assay == 'ELN'"

        ## remove FOSB, IL1A and IL1B as unreliable measurements
        pattern2 = "Assay in ['IL1A', 'IL1B', 'FOSB', 'EGR2', 'ELN']"

        #cav1 pattern needs finishing and the experiment object needs overwriting
        pattern3 = "cell_line in ['G', 'D', 'I'] and " \
                   "treatment == 'Control' and " \
                   "Assay == 'ELN'"

        pattern4 = "cell_line == 'I' and "\
                   "Assay == 'CAV1' and " \
                   "time == 8 and " \
                   "replicate == 3"


        df = remove(self.ddct, pattern1)
        df = remove(self.ddct, pattern2)
        df = remove(self.ddct, pattern3)
        df = remove(self.ddct, pattern4)

        return df


    def get_detailed_ddct_data(self, assay):
        """
        get data for just one gene (assay) to closer inspect data
        :param assay:
        :return:
        """
        df = self.ddct.query("Assay == '{}'".format(assay)).reset_index(drop=True)
        return pandas.pivot_table(df, values=0, columns='time', index=['cell_line', 'treatment', 'replicate'])

    def get_data_copasi_style_d3ct(self, gene, cell_line='A', plot=False, new_name=None,
                 filename=None, indep_dct=None):
        """
        Get fold change averaged data for gene. This is useful as
        it formats data for input into copasi
        :param exp:
        :param gene:
        :param cell_line:
        :param plot:
        :param new_name:
        :param filename:
        :return:
        """
        if gene not in self.genes:
            raise ValueError('{} is not in list of genes "{}"'.format(gene, self.genes))
        df = self.d3ct.query('cell_line == "{}" and Assay == "{}"'.format(cell_line, gene))
        df = pandas.pivot_table(df, values=0, index=['cell_line', 'replicate'], columns='time')
        mean = df.mean(axis=0)
        if new_name is not None:
            mean = pandas.DataFrame(mean, columns=[new_name])
        else:
            mean = pandas.DataFrame(mean, columns=[gene])

        mean = mean.reset_index()
        mean['time'] = mean['time']*60
        mean.loc[0] = [0, 0]

        if indep_dct is not None:
            for i in indep_dct:
                if i[-5:] != '_indep':
                    raise ValueError('Keys in indep_dct shuold have _indep as suffix')

                mean[i] = indep_dct[i]

        if filename is not None:
            mean.to_csv(filename, sep='\t', index=False)

        if plot:
            mean.plot(x='time')
            plt.show()

        return mean

    def get_data_copasi_style_dct(self, gene, cell_line='A', plot=False, new_name=None,
                                   filename=None, indep_dct=None):
        """
        Get fold change averaged data for gene. This is useful as
        it formats data for input into copasi
        :param exp:
        :param gene:
        :param cell_line:
        :param plot:
        :param new_name:
        :param filename:
        :return:
        """
        if gene not in self.genes:
            raise ValueError('{} is not in list of genes "{}"'.format(gene, self.genes))
        treat = self.treatment_data.query('cell_line == "{}" and Assay == "{}"'.format(cell_line, gene))['dct']
        treat.index = treat.index.droplevel('Assay')
        # treat = pandas.pivot_table(treat, index=['cell_line', 'treatment', 'replicate'], columns='time')
        treat = treat.groupby(level=['cell_line', 'treatment']).aggregate(numpy.mean).transpose()[cell_line]

        base = self.baseline_data.query('cell_line == "{}" and Assay == "{}"'.format(cell_line, gene))['dct']

        base = base.groupby(level=['cell_line', 'treatment']).aggregate(numpy.mean).transpose()[cell_line]
        baseline0 = pandas.DataFrame(base)['Baseline'][0]
        control = pandas.DataFrame(treat['Control'])
        tgf = pandas.DataFrame(treat['TGFb'])
        control.loc[0] = baseline0
        tgf.loc[0] = baseline0
        if new_name is not None:
            control.columns = [new_name]
            tgf.columns = [new_name]
            gene = new_name

        else:
            control.columns = [gene]
            tgf.columns = [gene]



        control = control.sort_index(level=0)
        tgf = tgf.sort_index(level=0)

        control = control.reset_index()
        tgf = tgf.reset_index()

        control['time'] = control['time'] * 60
        tgf['time'] = tgf['time'] * 60

        if indep_dct is not None:
            if len(indep_dct) != 2:
                raise ValueError('len of indep_dct should be 2 for this function. First key should be Control second TGFb')

            if 'Control' not in indep_dct.keys():
                raise ValueError

            if 'TGFb' not in indep_dct.keys():
                raise ValueError

            for i in indep_dct['Control']:
                if i[-6:] != '_indep':
                    raise ValueError('Keys in indep_dct shuold have _indep as suffix')

                if indep_dct['Control'][i] == '0_time_point_indep':
                    control[i] = float(control.loc[0, gene])

                else:
                    control[i] = indep_dct['Control'][i]

            for i in indep_dct['TGFb']:
                if i[-6:] != '_indep':
                    raise ValueError('Keys in indep_dct shuold have _indep as suffix')

                if indep_dct['TGFb'][i] == '0_time_point_indep':
                    tgf[i] = float(tgf.loc[0, gene])

                else:
                    tgf[i] = indep_dct['TGFb'][i]


                # tgf[i] = indep_dct['TGFb'][i]

        if filename is not None:
            ctrl_filename = "{}{}.txt".format(filename[:-4], '_control')
            tgf_filename = "{}{}.txt".format(filename[:-4], '_tgfb')
            control.to_csv(ctrl_filename, sep='\t', index=False)
            tgf.to_csv(tgf_filename, sep='\t', index=False)

        if plot:
            seaborn.set_context('talk', font_scale=2)
            seaborn.set_style('white')
            fig, ax = plt.subplots()
            plt.plot(tgf.time/60, tgf[gene], label='tgf', marker='o')
            plt.plot(control.time/60, control[gene], label='control', marker='o')
            plt.xlabel('Time(h)')
            plt.ylabel('dct')
            plt.legend(loc='best')
            seaborn.despine(fig=fig, top=True, right=True)
            plt.title('{}, {}'.format(g, cell_line))

            plt.show()

        return control, tgf

    @cached_property
    def design(self):
        df = pandas.read_csv(self._design)
        df = df.rename(columns={
            'Sub.experiment': 'sub_experiment',
            'Cell.ID': 'cell_id',
            'Cell.Line': 'cell_line',
            'Treatment': 'treatment',
            'Replicate': 'replicate',
            'Time.Point': 'time_point'
        })
        df = df.set_index(['sub_experiment', 'cell_id', 'treatment', 'replicate', 'time_point'])
        return df.sort_index(level=[0, 1, 2, 3, 4])

    def create_subexperiments(self):
        """
        Split design by subexperiment (1, 2, 3) and create SubExperiment
        objects
        :return:
        """
        subexperiments = {}
        for label, df in self.design.groupby(level=0):
            subexperiments[label] = SubExperiment(label, df.loc[label], self.root)
        return subexperiments

    @cached_property
    def treatment_data(self):
        """

        :return:
        """
        df = pandas.concat([i.treatment_data for i in list(self.subexperiments.values())])
        # df = df.swaplevel(0, 1)
        df = df.dropna(axis=0, how='all', thresh=self.max_nan)
        return df.sort_index(level=[0, 1, 2, 3])


    @cached_property
    def baseline_data(self):
        """

        :return:
        """
        df = pandas.concat([i.baseline_data for i in list(self.subexperiments.values())])
        # df = df.swaplevel(0, 1)
        return df.sort_index(level=[0, 1, 2, 3])

    def calculate_ddct(self, baseline_time=0):
        """
        For the time being this is TGF/baseline and control/baseline
        The choice of 0 or 96 will make a difference
        :return:
        """

        if baseline_time == 0:
            baseline = self.baseline_data['dct'][0]

        elif baseline_time == 96:
            baseline = self.baseline_data['dct'][96]

        elif baseline_time == 'average':
            baseline = (self.baseline_data['dct'][0] + self.baseline_data['dct'][96]) / 2

        control =self.treatment_data.query('treatment == "Control"')['dct']
        tgf = self.treatment_data.query('treatment == "TGFb"')['dct']
        control.index = control.index.droplevel(1)
        tgf.index = tgf.index.droplevel(1)
        baseline.index = baseline.index.droplevel(1)

        time = control.columns

        tgf = pandas.concat([tgf[i]/baseline for i in tgf.columns], axis=1)
        control = pandas.concat([control[i]/baseline for i in control.columns], axis=1)

        control.columns = time
        tgf.columns = time

        control = pandas.DataFrame(control.stack())
        tgf = pandas.DataFrame(tgf.stack())
        control['treatment'] = 'Control'
        tgf['treatment'] = 'TGFb'
        df = pandas.concat([control, tgf]).reset_index()
        df.rename(columns={0: 'ddct'})
        return df

    def calculate_d3ct(self):
        """
        calculate time course ddct divided by treated ddct
        :return:
        """
        data = deepcopy(self.ddct)
        data = data.set_index(['cell_line', 'replicate', 'Assay', 'time', 'treatment'])
        control = data.query('treatment == "Control"')#.reset_index(drop=True)
        tgfb = data.query('treatment == "TGFb"')#.reset_index(drop=True)
        control.index = control.index.droplevel(4)
        tgfb.index = tgfb.index.droplevel(4)
        return tgfb / control

    # def calculate_d4ct(self):
    #     """
    #     calculate time course ddct divided by treated ddct
    #     :return:
    #     """
    #     data = deepcopy(self.d3ct)
    #     data = data.reset_index()
    #     print data
    #     # data = data.set_index(['cell_line', 'replicate', 'Assay', 'time', 'treatment'])
    #     # control = data.query('treatment == "Control"')#.reset_index(drop=True)
    #     # tgfb = data.query('treatment == "TGFb"')#.reset_index(drop=True)
    #     # control.index = control.index.droplevel(4)
    #     # tgfb.index = tgfb.index.droplevel(4)
    #     # return tgfb / control

    def calculate_d4ct(self, pairs_dct=None):
        """
        calculate time course ddct divided by treated ddct
        :param pairs_dct:
            dict[cell_line] = Neonatal cell line
        :return:
        """
        #a, b, c, d, e, f, g, h, i
        if pairs_dct is None:
            pairs_dct={
                'D': 'A',
                'G': 'A',
                'E': 'B',
                'H': 'B',
                'F': 'C',
                'I': 'C',
            }
        data = deepcopy(self.d3ct)
        results_dct = {}
        for cell_line, normalizer in list(pairs_dct.items()):
            cell_data = data.query('cell_line =="{}"'.format(cell_line))
            norm_data = data.query('cell_line =="{}"'.format(normalizer))
            cell_data.index = cell_data.index.droplevel('cell_line')
            norm_data.index = norm_data.index.droplevel('cell_line')
            results_dct[r"$\frac{{{}}}{{{}}}$".format(cell_line, normalizer)] = cell_data / norm_data
        df = pandas.concat(results_dct)
        return df

class SubExperiment(object):
    """
    Class to hold a sub experiment formed of
    a HDFn, IR and A sub group.
    """
    def __init__(self, id, design, plate_directory):
        """
        :param id:
            Either 1, 2 or 3

        :param design:
            The submatrix containing the design for the
            sub-experiment

        :param plate_directory:
            location of folder containing the *well_data* files.
        """
        self._id = id
        self.design = design
        self.plate_directory = plate_directory

    def __str__(self):
        return "SubExperiment(id={}, cell_id={}, cell_lines='{}')".format(
            self.id, self.cell_id, self.cell_lines
        )

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, item):
        return self.plates[item]

    @cached_property
    def data(self):
        """
        get all the data, design plus raw plus normalized
        for subexperiment
        :return:
        """
        return pandas.concat([i.data for i in list(self.plates.values())])

    @cached_property
    def treatment_data(self):
        """

        :return:
        """
        return pandas.concat(i.treatment_data for i in list(self.plates.values()))

    @cached_property
    def baseline_data(self):
        """

        :return:
        """
        return pandas.concat(i.baseline_data for i in list(self.plates.values()))

    @cached_property
    def plates(self):
        plates = {}
        for label, df in self.design.groupby(by='Batch'):
            filename = os.path.join(self.plate_directory, df['Filename'].unique()[0])
            plates[label] = Plate(filename, df[df.Batch == label])
        return plates

    @cached_property
    def cell_id(self):
        return list(set(self.design.index.get_level_values(0)))

    @cached_property
    def treatments(self):
        return list(set(self.design.index.get_level_values(1)))

    @cached_property
    def replicates(self):
        return list(set(self.design.index.get_level_values(2)))

    @cached_property
    def time_points(self):
        return list(set(self.design.index.get_level_values(3)))

    @cached_property
    def id(self):
        return self._id

    @cached_property
    def cell_lines(self):
        return self.design['cell_line'].unique()

    @cached_property
    def lot_numbers(self):
        return self.design['Lot.Number'].unique()

    def get(self, query):
        """
        Get a subset of the design.

        A nice method but not actually used elsewhere.

        A query looks like this:
        query = {
            'treatment': ['Control', 'TGFb'],
            'cell_id': ['A', 'D', 'G'],
            'time_point': [0.5, 1, 2],
            'replicate': [1]
        }

        If any of these four elements are missing
        the default is to use all available values
        for that element.

        :param query:
            dictionary. keys: ['time_point', 'treatment',
                               'cell_id', 'replicate']
                        Values: lists of time points, treatments,
                                cell_ids or replicates
        :return:
            pandas.DataFrame containing subset of data
        """

        if not isinstance(query, dict):
            raise ValueError('Query should be a python dict')

        ##set default query containing all values
        default_query = {
            'cell_id': self.cell_id,
            'treatment': self.treatments,
            'time_point': self.time_points,
            'replicate': self.replicates
        }
        default_query.update(query)

        for k, v in list(query.items()):
            if not isinstance(v, list):
                query[k] = [v]


        available_labels = ['treatment', 'cell_id', 'time_point', 'replicate']


        labels = list(default_query.keys())
        if len(labels) not in list(range(1, 5)):
            raise ValueError('labels must be a list of length 1, 2, 3, or 4')

        for label in labels:
            if label not in available_labels:
                raise ValueError('"{}" not in accepted list of labels. These are '
                                 'accepted "{}"'.format(label, available_labels))


            if label is 'cell_id':
                cell_id_query = reduce(
                    lambda x, y: "{} or {}".format(x, y),
                    ['cell_id == "{}"'.format(i) for i in default_query[label]]
                )
            if label is 'treatment':
                treatment_query = reduce(
                    lambda x, y: "{} or {}".format(x, y),
                    ['treatment == "{}"'.format(i) for i in default_query[label]]
                )

            if label is 'time_point':
                time_query = reduce(
                    lambda x, y: "{} or {}".format(x, y),
                    ['time_point == {}'.format(i) for i in default_query[label]]
                )

            if label is 'replicate':
                replicate_query = reduce(
                    lambda x, y: "{} or {}".format(x, y),
                    ['replicate == {}'.format(i) for i in default_query[label]]
                )

        final_query = "({}) and ({}) and ({}) and ({})".format(
            cell_id_query, replicate_query, treatment_query,
            time_query
        )
        df = self.design.query(final_query)
        return self.data[self.data['Sample'].isin(df['Sample'])]



class Plate(object):
    """
    Class to represent a WaferGen plate. Each plate is 72 * 72 in
    dimension. Each column contains a specific gene while each row
    has a specific sample.
    """
    def __init__(self, filename, design):
        """
        :param filename:
            Name of data file containing the plate data

        :param design:
            Design corresponding to the plate.
        """
        self.filename = filename
        self.design = design
        self.id = os.path.split(filename)[1]
        self.columns = ['Row', 'Column', 'Assay',
                        'Sample', 'Ct']
        self.data = self._data()

        self.data = self.normalized_data

        self.baseline_data = self.organize_baseline()
        self.treatment_data = self.organize_treatments()


    def __str__(self):
        return 'Plate(batch={}, replicate={})'.format(
            self.batch, self.replicate
        )

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, item):
        return self.samples[item]

    def _data(self):
        """

        :return:
        """
        df = pandas.read_csv(self.filename, sep='\t')[self.columns]
        treatment, time, cell_line, rep = list(zip(*list(df['Sample'].str.split('_'))))
        df['treatment'] = treatment
        df['time'] = time
        df['cell_line'] = cell_line
        df['replicate'] = rep
        return df

    @cached_property
    def cell_id(self):
        return sorted(list(set(self.design.index.get_level_values(0))))

    @cached_property
    def treatments(self):
        return list(set(self.design.index.get_level_values(1)))

    @cached_property
    def replicate(self):
        return list(set(self.design.index.get_level_values(2)))[0]

    @cached_property
    def time_points(self):
        return list(set(self.design.index.get_level_values(3)))

    @cached_property
    def cell_line(self):
        return self.design['Cell.Line'].unique()

    @cached_property
    def lot_numbers(self):
        return self.design['Lot.Number'].unique()

    @cached_property
    def sub_experiment(self):
        return self.design['Sub.experiment'].unique()

    @cached_property
    def batch(self):
        return self.design['Batch'].unique()[0]

    @cached_property
    def samples(self):
        """
        Organize chip into readings per sample
        :return:
        """
        sample = {}
        for label, df in self.data.groupby(by='Sample'):
            sample[label] = Sample(label, df)
        return sample

    @cached_property
    def normalized_data(self):
        return pandas.concat([i.data for i in list(self.samples.values())])

    @cached_property
    def baseline(self):
        """
        Get baseline data without treatment time courses
        :return:
        """
        return self.data[self.data['treatment'] == 'Baseline']

    @cached_property
    def treatments(self):
        """
        Get treatment time courses without baseline
        data
        :return:
        """
        return self.data[self.data['treatment'] != 'Baseline']


    def organize_treatments(self):
        """

        :return:
        """
        self.treatments.loc[:, 'time'] = pandas.to_numeric(self.treatments['time'])
        self.treatments.loc[:, 'replicate'] = pandas.to_numeric(self.treatments['replicate'])
        df = self.treatments.pivot_table(
            index=['cell_line', 'treatment',
                   'replicate', 'Assay'],
            columns='time',
            values=['Sample', 'Ct', 'dct']
        )
        df['Ct'].columns = [float(i) for i in df['Ct'].columns]
        return df.sort_index(level=[0, 1], axis=1)

    def organize_baseline(self):
        """

        :return:
        """
        self.baseline.loc[:, 'time'] = pandas.to_numeric(self.baseline['time'])
        self.baseline.loc[:, 'replicate'] = pandas.to_numeric(self.baseline['replicate'])


        # pandas.to_numeric(self.baseline['time'])
        df = self.baseline.pivot_table(
            index=['cell_line', 'treatment',
                   'replicate', 'Assay'],
            columns='time',
            values=['Sample', 'Ct', 'dct']
        )
        df['Ct'].columns = [float(i) for i in df['Ct'].columns]
        return df.sort_index(level=[0, 1], axis=1)




class Sample(object):
    """
    Class to represent a single sample and data from all genes
    that belongs to that sample.

    Samples are created by the Plate class.

    """
    def __init__(self, id, data, reference_genes=['PPIA']):
        """

        :param id:
            Sample ID. Looks like cell_line_time_treatment_replicate

        :param data:
            Data belonging to this sample. This is passed down from
            the Plate class

        :param reference_genes:
            List of strings. Which reference genes to normalize to.
            The rest of the samples are normalized to the geometric
            mean of these genes.
        """
        self.id = id
        self.data = data
        self.reference_genes = reference_genes
        self.data = self.normalize


    def __str__(self):
        return "Sample(cell_line_type='{}', cell_line_name='{}', treatment='{}', " \
               "replicate={}, time={}h)".format(
            self.cell_line_type, self.cell_line_name,
            self.treatment, self.replicate, self.time,
        )

    def __repr__(self):
        return self.__str__()

    @cached_property
    def treatment(self):
        t = self.data['treatment'].unique()
        if not isinstance(t, (list, numpy.ndarray)):
            raise ValueError('Should be list. Got {}'.format(type(t)))

        if not len(t) == 1:
            raise ValueError('list should be length 1')

        return t[0]

    @cached_property
    def time(self):
        t = self.data['time'].unique()
        if not isinstance(t, (list, numpy.ndarray)):
            raise ValueError('Should be list. Got {}'.format(type(t)))

        if not len(t) == 1:
            raise ValueError('list should be length 1')

        return float(t[0])

    @cached_property
    def cell_line_type(self):
        t = self.data['cell_line'].unique()
        if not isinstance(t, (list, numpy.ndarray)):
            raise ValueError('Should be list. Got {}'.format(type(t)))

        if not len(t) == 1:
            raise ValueError('list should be length 1')

        return self.cell_line_lookup()[t[0]]

    @cached_property
    def cell_line_name(self):
        t = self.data['cell_line'].unique()
        if not isinstance(t, (list, numpy.ndarray)):
            raise ValueError('Should be list. Got {}'.format(type(t)))

        if not len(t) == 1:
            raise ValueError('list should be length 1')

        return t[0]

    @cached_property
    def replicate(self):
        t = self.data['replicate'].unique()
        if not isinstance(t, (list, numpy.ndarray)):
            raise ValueError('Should be list. Got {}'.format(type(t)))

        if not len(t) == 1:
            raise ValueError('list should be length 1')

        return int(t[0])

    @staticmethod
    def cell_line_lookup():
        return {'A': 'HDFn',
                'B': 'HDFn',
                'C': 'HDFn',
                'D': 'Senescent',
                'E': 'Senescent',
                'F': 'Senescent',
                'G': 'Adult',
                'H': 'Adult',
                'I': 'Adult',
                }

    @cached_property
    def genes(self):
        genes = self.data['Assay'].unique()
        return sorted(genes)

    @cached_property
    def normalize(self):
        df = pandas.concat([self.data[self.data['Assay'] == i] for i in self.reference_genes])
        geomean = gmean(df['Ct'])
        norm = pandas.DataFrame(2**(-(self.data['Ct'] - geomean)))
        norm.columns = ['dct']
        df = pandas.concat([self.data, norm], axis=1)
        return df.sort_values(by='Assay')



class Query(Experiment):
    """

    Query the experiment

    """
    def __init__(self, experiment, treatment, cell_id=None,
                 replicate=None, time=None,
                 gene=None, normed=False, averaged=False):
        """

        :param experiment:
            A Experiment obj

        :param treatment:
            Any of Baseline, TGFb or Control

        :param cell_id:
            Default: All.
            Any of A, B, C, D, E, F, G, H or I

        :param replicate:
            Default: all.
            Any of 1, 2, 3, 4, 5, 6

        :param time:
            if treatment is Baseline:
                0, 96
            else:
                0.5, 1, 2, 3, 4, 8, 12, 24, 48, 72, 96
            default: All

        :param gene:
            Any gene. Default all. Give incorrect gene
            and Python will present you with a list
            of valid genes


        :param normed:
            Whether to get normalized or raw data

        :param averaged:
            Not yet implemented
        """
        self.experiment = experiment
        self.cell_id = cell_id
        self.treatment = treatment
        self.replicate = replicate
        self.time = time
        self.gene = gene
        self.normed = normed
        self.averaged = averaged

        if not isinstance(self.experiment, Experiment):
            raise ValueError

        if self.treatment is 'Baseline':
            self.data = self.experiment.baseline_data
        else:
            self.data = self.experiment.treatment_data

        self.do_checks()

    def __str__(self):
        return "Query(\n\ttreatment='{}', normed='{}', \n\tcell_id='{}', " \
               "\n\treplicate='{}', \n\ttime='{}', \n\tgene='{}'\n)".format(
            self.treatment, self.normed, self.cell_id,
            self.replicate, self.time, self.gene
        )



    def do_checks(self):
        """
        verify integrity of user input
        :return:
        """
        # ## get valid experiment variables
        all_subexperiments = [1, 2, 3]
        all_plates = list(range(1, 19))
        all_cell_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
        all_samples = list(self.experiment.design['Sample'])
        all_genes = self.experiment.subexperiments[1].plates[1].samples[all_samples[0]].genes
        all_replicates = list(range(1, 7))
        all_time = [0.5, 1.0, 2.0, 3.0, 4.0, 8.0, 12.0, 24.0, 48.0, 72.0, 96.0]

        if self.time is None:
            if self.treatment is 'Baseline':
                self.time = [0.0, 96.0]
            else:
                self.time = all_time

        if self.cell_id is None:
            self.cell_id = all_cell_ids

        if self.gene is None:
            self.gene = all_genes

        if self.replicate is None:
            self.replicate = all_replicates

        if self.treatment is None:
            raise ValueError('treatment cannot be None. Specify one of "TGFb", "Control", "Baseline"')

        if not isinstance(self.treatment, str):
            raise ValueError('treatment must be a string. Got "{}" a "{}"'.format(
                self.treatment, type(self.treatment)
            ))

        if not isinstance(self.normed, bool):
            raise ValueError('normed argument should be boolean. Got "{}"'.format(
                type(self.normed)
            ))

        if not isinstance(self.time, list):
            self.time = [self.time]

        for time_point in self.time:
            if time_point not in sorted(list(set(self.data.columns.get_level_values(1)))):
                raise ValueError('"{}" is invalid time point. Valid time '
                                 'points are: {}'.format(
                    time_point, list(self.data.columns))
                )

    @property
    def result(self):
        if self.normed:
            return self.data.loc[self.treatment, self.cell_id,
                                 self.replicate, self.gene]['dct'][self.time]
        else:
            return self.data.loc[self.treatment, self.cell_id,
                                 self.replicate, self.gene]['Ct'][self.time]



if __name__ == '__main__':
    dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
    design_file = os.path.join(dire, 'new_design.csv')

    E = Experiment(design_file, outliers=False, )

    # eln = E.get_detailed_ddct_data('ELN')
    # print eln.to_csv('/home/b3053674/Documents/LargeStudy/SavedObjects/ELNdata.csv')
    # import pickle
    #
    experiment_pickle = os.path.join('/home/b3053674/Documents/LargeStudy/SavedObjects', 'Experiment.pickle')
    with open(experiment_pickle, 'wb') as f:
        pickle.dump(E, f)

    gene_list = ['ID1', 'SMAD7', 'COL1A1', 'ACTA2',
                 'COL1A2', 'CTGF', 'FN1', 'SERPINE1',
                 'JUNB']

    cell_lines = ['C', 'I']
    new_name = {
        'COL1A1': 'COL1A1mRNA',
        'COL1A2': 'COL1A2mRNA',
        'SMAD7': 'Smad7mRNA',
        'CTGF': 'CTGFmRNA',
    }
    for g in gene_list:
        for c in cell_lines:
            fname = os.path.join('/home/b3053674/Documents/Models/2018/05_May/TGFModel/Data', '{}_{}.txt'.format(c, g))
            if g in new_name.keys():
                new = new_name[g]

            else:
                new = g

            fig = E.get_data_copasi_style_dct(
                g, cell_line=c, plot=True,
                filename=fname,
                indep_dct={'TGFb':        {'TGFb_indep': 1,
                                           '{}_indep'.format(new): '0_time_point_indep',
                                           'Neonatal_indep': 1 if c == 'C' else 0},
                           'Control':     {'TGFb_indep': 0,
                                           '{}_indep'.format(new): '0_time_point_indep',
                                           'Neonatal_indep': 1 if c == 'C' else 0},
                           },
                new_name=new
            )


    # g2 = E.get_data_copasi_style_dct('COL1A1', cell_line='G', plot=False, new_name=None, filename=None)
    # print(g)
    # print(g2)
    # print(E.baseline_data.to_csv(r'/home/b3053674/Documents/LargeStudy/timeseriesAnalysis/baseline_data.csv'))

    # print(E.d3ct.to_csv('/home/b3053674/Documents/LargeStudy/timeseriesAnalysis/'))
    #
    # with open(experiment_pickle, 'r') as f:
    #     E = pickle.load(f)

    # E.ddct.to_csv(r'/home/b3053674/Documents/LargeStudy/SavedObjects/ddct_data10-03-2018.csv', index=False)
    # E.d3ct.reset_index().to_csv(r'/home/b3053674/Documents/LargeStudy/SavedObjects/d3ct_data10-03-2018.csv', index=False)

    # print E.treatment_data



    # df = pandas.concat([E.treatment_data.stack(), E.baseline_data.stack()])
    # df.to_csv(os.path.join(dire, 'DataFromWaferGen2.csv'))
    # col1a1 = E.get_gene('COL1A1', plot=False, new_name='COL1A1mRNAObs',
    #                   filename=r'/home/b3053674/Documents/Models/2018/02_Feb/SmadErkModel/COL1A1mRNAObs.csv')
    # col1a2 = E.get_gene('COL1A2', plot=False, new_name='COL1A2mRNAObs',
    #                   filename=r'/home/b3053674/Documents/Models/2018/02_Feb/SmadErkModel/COL1A2mRNAObs.csv')
    # smad7 = E.get_gene('SMAD7', plot=False, new_name='Smad7mRNAObs',
    #                  filename=r'/home/b3053674/Documents/Models/2018/02_Feb/SmadErkModel/Smad7mRNAObs.csv')
    # ctgf = E.get_gene('CTGF', plot=False, new_name='CTGFmRNAObs',
    #                 filename=r'/home/b3053674/Documents/Models/2018/02_Feb/SmadErkModel/CTGFmRNAObs.csv')
    # pai1 = E.get_gene('SERPINE1', plot=False, new_name='PAI1mRNAObs',
    #                 filename=r'/home/b3053674/Documents/Models/2018/02_Feb/SmadErkModel/PAI1mRNAObs.csv')
    # id1 = E.get_gene('ID1', plot=False, new_name='ID1mRNAObs',
    #                filename=r'/home/b3053674/Documents/Models/2018/02_Feb/SmadErkModel/ID1mRNAObs.csv')
    # id1 = E.get_gene('ID1', plot=False, new_name='ID1mRNAObs',
    #                filename=r'/home/b3053674/Documents/Models/2018/02_Feb/SmadErkModel/ID1mRNAObs.csv')





























