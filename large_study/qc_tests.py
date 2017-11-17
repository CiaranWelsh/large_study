import os, pandas, numpy, glob, unittest
from parse import *
from qc import *





class TestQC(unittest.TestCase):
    def setUp(self):
        self.dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
        self.design_file = os.path.join(dire, 'new_design.csv')


    # def test_gene_batches(self):
    #     """
    #
    #     :return:
    #     """
    #
    #     gene_batches = get_gene_batches()
    #     self.assertEqual(len(gene_batches[0]) * len(gene_batches), 72)
    #
    #
    # def test_pca(self):
    #     exp = Experiment(self.design_file)
    #     for batch in get_gene_batches()[:2]:
    #         pca(exp, batch)

        # p = pca(exp)



if __name__ == '__main__':
    unittest.main()