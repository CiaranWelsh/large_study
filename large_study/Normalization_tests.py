import pandas
import numpy
import unittest
from . import parse
import os, glob


class TestQuantileNormalization(unittest.TestCase):
    """
    Using the example from the wikipedia
    page on normalization
    :return:
    """
    def setUp(self):
        self.df = pandas.DataFrame({'C1': {'A': 5, 'B': 2, 'C': 3, 'D': 4},
                               'C2': {'A': 4, 'B': 1, 'C': 4, 'D': 2},
                               'C3': {'A': 3, 'B': 4, 'C': 6, 'D': 8}})
        self.data_dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
        self.new_dire = r'/home/b3053674/Documents/LargeStudy/QuantileNormalizedData'

    def test_normalized_get_data_mean(self):
        """
        test that new quantile normalized data exists in data
        :return:
        """
        normed_data = parse.QuantileNormalization(self.data_dire, self.new_dire).normalized_data
        for i in normed_data:
            self.assertIn('CtQuantileNormalizedMean', normed_data[i].columns)

    def test_normalized_get_data_median(self):
        """
        test that new quantile normalized data exists in data
        :return:
        """
        normed_data = parse.QuantileNormalization(self.data_dire, self.new_dire).normalized_data
        for i in normed_data:
            self.assertIn('CtQuantileNormalizedMedian', normed_data[i].columns)

    def test_frame(self):
        """
        Answer:
            A    5.67    4.67    2.00
            B    2.00    2.00    3.00
            C    3.00    4.67    4.67
            D    4.67    3.00    5.67
        :return:
        """
        answer = pandas.DataFrame({'C1': {'A': 5.67, 'B': 2.00, 'C': 3.00, 'D': 4.67},
                                   'C2': {'A': 4.67, 'B': 2.00, 'C': 4.67, 'D': 3.00},
                                   'C3': {'A': 2.00, 'B': 3.00, 'C': 4.67, 'D': 5.67}})

        actual = parse.QuantileNormalization(self.data_dire, self.new_dire).normalize(self.df).as_matrix()
        numpy.testing.assert_almost_equal(actual, answer.as_matrix(), decimal=2)

    def test_to_file(self):
        parse.QuantileNormalization(self.data_dire, self.new_dire).to_file()
        self.assertTrue(os.path.isdir(self.new_dire))

    # def test_t(self):
    #     Q = parse.QuantileNormalization(self.data_dire, self.new_dire)
    #     # Q.normalize_data()
    #     # Q.to_file()
    #     print Q.plot_distributions('Ct')
    #     print Q.plot_distributions('QuantileNormalizedMedian')
    #     print Q.plot_distributions('QuantileNormalizedMean')
    #     # print Q.read_data2()
    #     # Q.qnorm_data()




if __name__ == '__main__':
    unittest.main()