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
import os, glob
from copy import deepcopy
from scipy.stats.mstats import gmean
import logging
from cached_property import cached_property

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
    def __init__(self, design, max_nan=6):
        """

        :param design:
            Path to the design file.

        :param max_nan:
            If max_nan NaN's present in treated or control profile,
            remove the data entry. default = 6.
        """
        self._design = design
        self.max_nan = max_nan
        self.root = os.path.dirname(design)
        self.subexperiments = self.create_subexperiments()
        ## nested unzipping of nested list comprehension to retrieve all sample names from all plates
        self.all_samples = reduce(
            lambda x, y: x + y,
            [self.subexperiments[i].plates[j].samples.keys() for i in self.subexperiments for j in self.subexperiments[i].plates]
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

    def __getitem__(self, item):
        return self.subexperiments[item]

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
        df = pandas.concat([i.treatment_data for i in self.subexperiments.values()])
        # df = df.swaplevel(0, 1)
        df = df.dropna(axis=0, how='all', thresh=self.max_nan)
        return df.sort_index(level=[0, 1, 2, 3])


    @cached_property
    def baseline_data(self):
        """

        :return:
        """
        df = pandas.concat([i.baseline_data for i in self.subexperiments.values()])
        # df = df.swaplevel(0, 1)
        return df.sort_index(level=[0, 1, 2, 3])

#
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
        return pandas.concat([i.data for i in self.plates.values()])

    @cached_property
    def treatment_data(self):
        """

        :return:
        """
        return pandas.concat(i.treatment_data for i in self.plates.values())

    @cached_property
    def baseline_data(self):
        """

        :return:
        """
        return pandas.concat(i.baseline_data for i in self.plates.values())

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

        for k, v in query.items():
            if not isinstance(v, list):
                query[k] = [v]


        available_labels = ['treatment', 'cell_id', 'time_point', 'replicate']


        labels = default_query.keys()
        if len(labels) not in range(1, 5):
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
        treatment, time, cell_line, rep = zip(*list(df['Sample'].str.split('_')))
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
        return pandas.concat([i.data for i in self.samples.values()])

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
            values=['Sample', 'Ct', 'Norm2ref']
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
            values=['Sample', 'Ct', 'Norm2ref']
        )
        df['Ct'].columns = [float(i) for i in df['Ct'].columns]
        return df.sort_index(level=[0, 1], axis=1)




class Sample(object):
    """
    Class to represent a single sample and data from all genes
    that belongs to that sample.

    Samples are created by the Plate class.

    """
    def __init__(self, id, data, reference_genes=['B2M', 'PPIA']):
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
        norm.columns = ['Norm2ref']
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
        all_plates = range(1, 19)
        all_cell_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
        all_samples = list(self.experiment.design['Sample'])
        all_genes = self.experiment.subexperiments[1].plates[1].samples[all_samples[0]].genes
        all_replicates = range(1, 7)
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

        if not isinstance(self.treatment, (str, unicode)):
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
                                 self.replicate, self.gene]['Norm2ref'][self.time]
        else:
            return self.data.loc[self.treatment, self.cell_id,
                                 self.replicate, self.gene]['Ct'][self.time]



if __name__ == '__main__':
    dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
    design_file = os.path.join(dire, 'new_design.csv')

    E = Experiment(design_file)
    print E.treatment_data
    # res = Query(E, 'Control', cell_id='A').result
    # print res