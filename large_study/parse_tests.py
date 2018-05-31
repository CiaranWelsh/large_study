import os, glob
import pandas, numpy
import unittest
from .parse import *



class TestExperiment(unittest.TestCase):
    def setUp(self):
        self.design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'
        self.design_file = r'/home/b3053674/Documents/LargeStudy/QuantileNormalizedData/new_design.csv'
        self.E = Experiment(self.design_file)

    def test_number_of_samples(self):
        """
        should be 1296
        :return:
        """
        self.assertEqual(1296, len(self.E.all_samples))

    def test_subexperiments1(self):
        """
        Test that three sub experiments
        have been created
        :return:
        """
        self.assertEqual(len(self.E.subexperiments), 3)

    def test_subexperiments2(self):
        """
        Test that three sub experiments
        have been created with correct id's
        :return:
        """
        self.assertEqual([1, 2, 3], [i.id for i in list(self.E.subexperiments.values())])

    def test_design_length(self):
        """
        Make sure design file accounts for all
        samples
        :return:
        """
        self.assertEqual(self.E.design.shape[0], 1296)

    def test_design_keys(self):
        """
        Make sure design file accounts for all
        samples
        :return:
        """
        l = ['cell_line', 'Lot.Number',
             'Batch', 'Filename', 'WG.Plate',
             'Sample', 'Treatment Start Date']
        self.assertTrue([i in list(self.E.design.keys()) for i in l])

    def test_treatment_data(self):
        cell_lines = sorted(list(set(self.E.treatment_data.index.get_level_values('cell_line'))))
        self.assertListEqual(cell_lines, ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'])

    def test_baseline_data(self):
        cell_lines = sorted(list(set(self.E.baseline_data.index.get_level_values('cell_line'))))
        self.assertListEqual(cell_lines, ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'])

    def test_max_nan(self):
        exp = Experiment(self.design_file, max_nan=12)
        self.assertEqual(exp.treatment_data.shape, (7381, 22))

    def test_max_nan2(self):
        exp = Experiment(self.design_file, max_nan=8)
        self.assertEqual(exp.treatment_data.shape, (7302, 22))

    def tearDown(self):
        del self.E


class TestSubExperiment(unittest.TestCase):
    def setUp(self):
        large_study_dir = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'

        self.design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'
        self.design_file = r'/home/b3053674/Documents/LargeStudy/QuantileNormalizedData/new_design.csv'

        self.E = Experiment(self.design_file)
        self.sub1 = self.E.subexperiments[1]

    def test_sub1_cell_id(self):
        """
        Test contents of sub experiment 1
        :return:
        """
        cell_ids = ['A', 'D', 'G']
        self.assertEqual(self.sub1.cell_id, cell_ids)

    def test_sub1_treatments(self):
        """
        Test contents of sub experiment 1
        :return:
        """
        treatments = ['Control', 'TGFb', 'Baseline']
        self.assertEqual(self.sub1.treatments, treatments)


    def test_plates_length(self):
        """
        6 plates per sub experiment
        :return:
        """
        self.assertEqual(len(self.sub1.plates), 6)

    def test_plates(self):
        """
        tests plates are created
        :return:
        """
        self.assertTrue([type(i) == Plate for i in self.sub1.plates])

    def test_well_data_sub_exp1(self):
        """
        make sure that sub experiment 1 has
        the correct cell lines in
        :return:
        """
        self.assertListEqual(
            self.E.subexperiments[1].plates[1].cell_id, ['A', 'D', 'G']
        )

    def test_well_data_sub_exp2(self):
        """
        make sure that sub experiment 2 has
        the correct cell lines in
        :return:
        """
        self.assertListEqual(
            self.E.subexperiments[2].plates[8].cell_id, ['B', 'E', 'H']
        )

    def test_well_data_sub_exp3(self):
        """
        make sure that sub experiment 3 has
        the correct cell lines in
        :return:
        """
        self.assertListEqual(
            self.E.subexperiments[3].plates[14].cell_id, ['C', 'F', 'I']
        )

    def test_data_cell_lines(self):
        cell_lines = self.sub1.data['cell_line'].unique()
        self.assertListEqual(sorted(cell_lines), ['A', 'D', 'G'])

    def test_treatment_data(self):
        cell_lines = list(set(self.sub1.treatment_data.index.get_level_values('cell_line')))
        self.assertListEqual(sorted(cell_lines), ['A', 'D', 'G'])

    def test_baseline_data(self):
        cell_lines = list(set(self.sub1.baseline_data.index.get_level_values('cell_line')))
        self.assertListEqual(sorted(cell_lines), ['A', 'D', 'G'])

    def test_all_data(self):
        """
        6 plates per sub experiment
        :return:
        """
        self.assertEqual(len(self.sub1.plates), 6)


class PlateTests(unittest.TestCase):
    def setUp(self):
        self.design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'
        self.design_file = r'/home/b3053674/Documents/LargeStudy/QuantileNormalizedData/new_design.csv'
        self.E = Experiment(self.design_file)
        self.sub1 = self.E.subexperiments[1]
        self.plate1 = self.sub1.plates[1]
        self.plate2 = self.sub1.plates[2]
        self.plate3 = self.sub1.plates[3]
        self.plate4 = self.sub1.plates[4]
        self.plate5 = self.sub1.plates[5]
        self.plate6 = self.sub1.plates[6]

    def test_plate1_filename(self):
        filename = '/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/GSS2375_WellData_p1.txt'
        self.assertEqual(os.path.split(self.plate1.filename)[1], os.path.split(filename)[1])

    def test_plate1_replicate(self):
        self.assertEqual(self.plate1.replicate, 1)

    def test_plate_well_data_repeat1(self):
        """
        test that plate1 contains only repeat 1 data
        :return:
        """
        replicate = '1'
        self.assertEqual(self.plate1.data.replicate.unique()[0], replicate)


    def test_plate_well_data_repeat2(self):
        """
        test that plate1 contains only repeat 1 data
        :return:
        """
        self.assertEqual(self.plate2.data.replicate.unique()[0], '2')

    def test_plate_well_data_repeat3(self):
        """
        test that plate1 contains only repeat 1 data
        :return:
        """
        self.assertEqual(self.plate3.data.replicate.unique()[0], '3')

    def test_plate_well_data_repeat4(self):
        """
        test that plate1 contains only repeat 1 data
        :return:
        """
        self.assertEqual(self.plate4.data.replicate.unique()[0], '4')

    def test_plate_well_data_repeat5(self):
        """
        test that plate1 contains only repeat 1 data
        :return:
        """
        self.assertEqual(self.plate5.data.replicate.unique()[0], '5')

    def test_plate_well_data_repeat6(self):
        """
        test that plate1 contains only repeat 1 data
        :return:
        """
        self.assertEqual(self.plate6.data.replicate.unique()[0], '6')

    def test_samples_length(self):
        self.assertEqual(len(self.plate1.samples), 72)

    def test_baseline(self):
        """
        Test that baseline can be separated from
        the control and treated
        :return:
        """
        treatments = list(set(self.plate1.baseline_data.index.get_level_values(level='treatment')))
        self.assertListEqual(treatments, ['Baseline'])
    #
    def test_treatments(self):
        """
        Test that baseline can be separated from
        the control and treated
        :return:
        """
        treatments = list(set(self.plate1.treatment_data.index.get_level_values(level='treatment')))
        self.assertListEqual(treatments, ['Control', 'TGFb'])


    def tearDown(self):
        pass


class TestSamples(unittest.TestCase):
    def setUp(self):
        self.design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'
        self.design_file = r'/home/b3053674/Documents/LargeStudy/QuantileNormalizedData/new_design.csv'
        self.E = Experiment(self.design_file)
        self.sub1 = self.E.subexperiments[1]
        self.plate1 = self.sub1.plates[1]
        self.s1 = self.plate1.samples['TGFb_24_A_1']

    def test_sample_normalization(self):
        """
        Compare manual normalization with
        that done by the Sample class.
        Normalization is:
            2**-(gene_i - reference)

        where:
            gene_i is a samples measurement
            of each each

            reference is the geometric mean
            of the reference genes (B2M and PPIA)
            in a sample
        :return:
        """
        B2M_sample1_plate1_sub_exp1 = 20.57081
        PPIA_sample1_plate1_sub_exp1 = 22.73551
        ACTA2_normed = 0.45694803996582645

        print(self.s1.data)
        self.assertAlmostEqual(float(self.s1.data.query('Assay == "ACTA2"')['Norm2ref']), ACTA2_normed)











if __name__ == '__main__':
    unittest.main()





