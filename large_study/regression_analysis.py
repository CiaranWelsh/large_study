from . import parse
import os, glob, pandas, numpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import seaborn
from scipy.stats.mstats import pearsonr

"""
Time Lag Regression
======================

I've been switching between Python and R but R is annoying so I've gone back to Python. Initially 
I created a couple of R functions to perform linear regression between a single predictor 
and a response variable (i.e. COL1A1 and CTGF). Each data set was split by treatments and cell
types in order to compare coefficients and pearsons correlation. This has worked fairly well. 

I now want to introduce a time lag parameter. The idea is that I want to correlate a predictor
variable at some time in the past, specified by dt with the response variable. 


Plan
====
1) isolate the two variables of interest (COL1A1 and CTGF for instance)
Do control, treated and treated/control. Not baseline. Therefore
remove baseline and added TGF/Control column

2) Interpolate the data using cubic spline for example. specify resolution.
3) add time delay to the predictor variable
4) Calculate pearsons coef. Also build in mutual information options for comparison
5) Find a way of maximizing pearsons or MI accross range of time delay


"""




class TimeLagRegression(object):
    def __init__(self, data, x, y, interpolation_kind='cubic', step_size=0.1):
        self.data = data
        self.x = x
        self.y = y
        self.interpolation_kind = interpolation_kind
        self.step_size = step_size

        ## remove any anomolies (id'ed with PCA)
        self.data = self.remove_anomolies()

        ## isolate the x and y variabels
        self.data = self.get_variables()

        ## remove baseline
        self.data = self.drop_baseline()

        ## unstack the data so that time is along columns
        self.data = self.unstack_data()

        ## calculate TGFb / control and add to dataframe
        self.data = self.calc_fc()

        ## interpolate the data
        self.interp_data = self.interpolate()

        self.do()

        # self.time_delay_data = self.calculate_time_delay_data(1)

    def remove_anomolies(self):
        """
        remove any anomolies such as repeat 6, TGFb, 8h cell line F
        :return:
        """
        return self.data.query("cell_line != 'F' or treatment != 'TGFb' or  replicate != 6 or time == 8")

    def get_variables(self):
        return self.data.query("Assay in ['{}', '{}']".format(self.x, self.y))

    def drop_baseline(self):
        """
        baseline has too few time points to interpolate and is
        thereby dropped.
        :return:
        """
        return self.data.query('treatment != "Baseline"')

    def unstack_data(self):
        """
        get data so that time is along columns
        :return:
        """
        df = self.data.unstack()
        return df['Norm2ref']

    def calc_fc(self):
        """
        into TGFb, control and TGFb/control
        :return:
        """
        data = self.data
        data.index = data.index.swaplevel(0, 1)
        fc = data.loc['TGFb'] / data.loc['Control']
        fc['treatment'] = ['TGFb / Control'] * fc.shape[0]
        fc = fc.reset_index()
        fc = fc.set_index(['treatment', 'cell_line', 'replicate', 'Assay'])
        return pandas.concat([data, fc])

    def interpolate1timecourse(self, x, y, **kwargs):
        """
        Interpolate a time course using scipy.interpolation.interp1d
        :param data: vector of data points to interpolate
        :param kwargs: passed on to interp1d
        :return:
        """
        f = interp1d(x, y, kind=kwargs.get('kind'))
        x_new = numpy.arange(x[0], x[-1], step=0.1)
        y_new = f(x_new)
        return pandas.DataFrame(y_new, index=x_new).transpose()

    def interpolate(self):
        df_list = []
        for label, df in self.data.groupby(level=['treatment', 'cell_line', 'replicate', 'Assay']):
            x = list(df.loc[label].index)
            y = df.loc[label].values
            interp = self.interpolate1timecourse(x, y, kind=self.interpolation_kind)
            interp.index = df.index
            interp.columns.name = 'time'
            df_list.append(interp)
        return pandas.concat(df_list)

    def plot_interpolated(self):
        """

        :return:
        """
        data = self.interp_data.stack().reset_index()
        print(data.head())
        x = data.query("Assay == '{}'".format(self.x))
        y = data.query("Assay == '{}'".format(self.y))
        print(x.shape)
        print(y.shape)
        for label, df, in data.groupby(by=['treatment', 'cell_line']):
            plt.figure()
            seaborn.tsplot(data=df, time='time', unit='replicate', condition='Assay', value=0)
            plt.title('{}_{}'.format(label[0], label[1]))
        plt.show()

    def calculate_time_delay_data(self, dt):
        """
        Example Data:

                time: 0, 1, 2, 3, 4h
                x:    4, 5, 6, 7, 8
                y:    7, 5, 4, 3, 2

                time: 0, 1, 2, 3, 4h
                x:    4, 5, 6, 7, 8
                y:    5, 4, 3, 2

        :param dt:
        :return:
        """
        x = self.interp_data.query('Assay == "{}"'.format(self.x))
        y = self.interp_data.query('Assay == "{}"'.format(self.y))

        ##introduce the time delay
        new_x = [i + dt for i in list(x.columns)]
        x.columns = new_x

        ## find intersection between x.columns and y.columns
        x_intersect_y = list(set(x.columns).intersection(set(y.columns)))

        x_delay = x[x_intersect_y]
        y = y[x_intersect_y]

        return x, y


    def do(self):
        """
        iterate over cell line and treatment.
        For each:
            calculate pearsons correlation
        :return:
        """

        x, y = self.calculate_time_delay_data(1)
        # print pearsonr(x, y)
        for label, df in x.groupby(level=['treatment', 'cell_line']):
            print(df)


    def pearsons_correlation(self, x, y):
        return pearsonr(x, y)








if __name__ == '__main__':
    dire = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant'
    dire = r'C:\Users\Ciaran\Documents\large_study\GSS2375_WB_NewDur_Grant'
    design_file = os.path.join(dire, 'new_design.csv')

    data_file = os.path.join(dire, 'DataFromWaferGen2.csv')

    df = pandas.read_csv(data_file, index_col=[0, 1, 2, 3, 4])
    df = pandas.DataFrame(df['Norm2ref'])

    TLR = TimeLagRegression(df, 'CTGF', 'COL1A1')


















