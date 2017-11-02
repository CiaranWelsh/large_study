import os
import pandas
import unittest
import plot
import parse

class QueryDataTests(unittest.TestCase):
    def setUp(self):
        self.design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'
        self.exp = parse.Experiment(self.design_file)
        # query = {
        #     'treatment': 'TGFb',
        #     'cell_id': ['A', 'D'],
        #     'subexperiment': 1,
        #     'replicate': [1,2,3],
        #     'gene': 'SMAD7',
        #     'time': 'all'
        # }


    def test_query_0(self):
        """
        Answer:
            Baseline  A         1         ACTA2      0.120495   0.212906
        :return:
        """
        Q = plot.Query(self.exp,
                       treatment='Baseline', normed=True)
        # self.assertEqual(Q.result[0])
        ans0 = float(Q.result.loc['Baseline', 'A', 1, 'ACTA2'][0])
        self.assertAlmostEqual(ans0, 0.12049530069298563)

    def test_query_96(self):
        """
        Answer:
            Baseline  A         1         ACTA2      0.120495   0.212906
        :return:
        """
        Q = plot.Query(self.exp,
                       treatment='Baseline', normed=True)
        # self.assertEqual(Q.result[0])
        ans96 = float(Q.result.loc['Baseline', 'A', 1, 'ACTA2'][96])
        self.assertAlmostEqual(ans96, 0.21290572556238796)


    def test(self):
        """
        Get
            * Gene: Smad7
            * Raw data
            * treatment: TGFb
            * time: 8h
            * cell line: C
            * repeat: 1

        Answer:
            ct: 25.60106
        :return:
        """
        Q = plot.Query(self.exp, treatment='TGFb', time=8,
                           cell_id='C', gene='SMAD7', replicate=1, normed=False)

        self.assertAlmostEqual(float(Q.result), 25.60106)










if __name__ == '__main__':
    unittest.main()
