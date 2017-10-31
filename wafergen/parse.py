import pandas, numpy
import os, glob
from copy import deepcopy
from scipy.stats.mstats import gmean
import logging
from cached_property import cached_property

# FORMAT = "%(asctime)s"
logging.basicConfig(level=logging.DEBUG)




class Sample(object):
    def __init__(self, id, data, reference_genes=['B2M', 'PPIA']):
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


class Plate(object):
    """
    Each plate is 72 * 72 in dimensions
    Each column contains a specific gene
    Each row has a specific sample
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
        # self.reference_genes = ['BGN', '36B4', 'PPIA']
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
        sample = []
        for label, df in self.data.groupby(by='Sample'):
            sample.append(Sample(label, df))
        return sample

    @cached_property
    def normalized_data(self):
        return pandas.concat([i.data for i in self.samples])

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
        return self.treatments.pivot_table(
            index=['cell_line', 'treatment',
                   'replicate', 'Assay'],
            columns='time',
            values=['Sample', 'Ct', 'Norm2ref']
        )

    def organize_baseline(self):
        """

        :return:
        """
        return self.baseline.pivot_table(
            index=['cell_line', 'treatment',
                   'replicate', 'Assay'],
            columns='time',
            values=['Sample', 'Ct', 'Norm2ref']
        )

class SubExperiment(object):
    """
    """
    def __init__(self, id, design, plate_directory):
        self._id = id
        self.design = design
        self.plate_directory = plate_directory

    def __str__(self):
        return "SubExperiment(id={}, cell_id={}, cell_lines='{}')".format(
            self.id, self.cell_id, self.cell_lines
        )

    def __repr__(self):
        return self.__str__()

    @cached_property
    def data(self):
        """
        get all the data, design plus raw plus normalized
        for subexperiment
        :return:
        """
        return pandas.concat([i.data for i in self.plates])

    @cached_property
    def treatment_data(self):
        """

        :return:
        """
        return pandas.concat(i.treatment_data for i in self.plates)

    @cached_property
    def baseline_data(self):
        """

        :return:
        """
        return pandas.concat(i.baseline_data for i in self.plates)

    @cached_property
    def plates(self):
        plates = []
        for label, df in self.design.groupby(by='Batch'):
            filename = os.path.join(self.plate_directory, df['Filename'].unique()[0])
            plates.append(Plate(filename, df[df.Batch == label]))
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
        Get a subset of the design

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
        return self.all_data[self.all_data['Sample'].isin(df['Sample'])]
        # print df[df['Sample'].isin(self.design['Sample'])]



class Experiment(object):
    """

    """
    def __init__(self, design):
        self._design = design
        self.root = os.path.dirname(design)
        self.subexperiments = self.create_subexperiments()

    @property
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
        subexperiments = []
        for label, df in self.design.groupby(level=0):
            subexperiments.append(SubExperiment(label, df.loc[label], self.root))
        return subexperiments

    @cached_property
    def treatment_data(self):
        """

        :return:
        """
        return pandas.concat([i.treatment_data for i in self.subexperiments])

    @cached_property
    def baseline_data(self):
        """

        :return:
        """
        return pandas.concat([i.baseline_data for i in self.subexperiments])

#
# if __name__ == '__main__':
#     large_study_dir = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
#
#
#
#
#     design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'
#     E = Experiment(design_file)
#     s1 = E.subexperiments[0]
#
#     query = {
#         'cell_id': ['A', 'D', 'e'],
#         'replicate': [1, 3, 5]
#     }
#     #(cell_line == "A" or cell_line == "D" or cell_line == "G") and (replicate == 1) and (treatment == "Control" or treatment == "TGFb") and (time_point == 0.5 or time_point == 1 or time_point == 2)
#
#     f = r'/home/b3053674/Documents/LargeStudy/df.csv'
#     p1 = s1.plates[3]#[5].data.to_csv(f)
    # print p1.data.replicate.unique()
    # print s1.plates
    # s1.all_data.to_csv(f)
    # print s1.get(query).to_csv(f)
    # print s1.get_data()
    # print s1.design.query('(cell_id == "A" or cell_id == "D" or cell_id == "G")')











    # dire = r'C:\Users\Ciaran\Documents\LargeStudy\GSS2375_WB_NewDur_Grant'
    # files = glob.glob(os.path.join(dire, '*WellData*'))
    # p = Plate(files[0])
    # print p.data


    # import copy
    # df = copy.deepcopy(E.design)

    # treat = list(df['Treatment'])
    # time = list(df['Time.Point'])
    # id = list(df['Cell.ID'])
    # rep = list(df['Replicate'])
    # sample = []
    # for i in range(len(treat)):
    #     if time[i] == 0.5:
    #         sample.append("{}_{}_{}_{}".format(
    #             treat[i], time[i], id[i], rep[i]
    #         )
    #     )
    #     else:
    #         sample.append("{}_{}_{}_{}".format(
    #             treat[i], int(time[i]), id[i], rep[i]
    #         ))
    # df['Sample'] = sample
    # print df.to_csv(r'/home/b3053674/Documents/LargeStudy/new_design.csv')
    # s1 = E.subexperiments[1]
    # print s1.design.head()
    # plate = s1.plates[0]
    # genes = plate.samples[0].genes
    # print plate

