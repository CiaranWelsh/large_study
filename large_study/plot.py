# -*- coding: utf8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from parse import *
import functools32


app = dash.Dash()

subexperiment = 1
plate = 1
cell_id = 'A'
treatments = 'TGFb'
time = 'all'
repeats = 'all'
averaged = False
normalized = False

design_file = r'/home/b3053674/Documents/LargeStudy/GSS2375_WB_NewDur_Grant/new_design.csv'
exp = Experiment(design_file)
sample_keys = exp.subexperiments[1].plates[1].samples.keys()
genes = exp.subexperiments[1].plates[1].samples[sample_keys[0]].genes


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





app.layout = html.Div(
    style={'textAlign': 'center'},
    children=[
        html.Div([
            dcc.Checklist(
                options=[
                    {'label': 'TGFb', 'value': 'TGFb'},
                    {'label': 'Control', 'value': 'Control'},
                    {'label': 'Baseline', 'value': 'Baseline'},
                ],
                id='treatment_checklist',
                values=['TGFb', 'Control', 'Baseline'],
            ),
        ]),

        dcc.Checklist(
            options=[
                {'label': 'A', 'value': 'A'},
                {'label': 'B', 'value': 'B'},
                {'label': 'C', 'value': 'C'},
                {'label': 'D', 'value': 'D'},
                {'label': 'E', 'value': 'E'},
                {'label': 'F', 'value': 'F'},
                {'label': 'G', 'value': 'G'},
                {'label': 'H', 'value': 'H'},
                {'label': 'I', 'value': 'I'},
            ],
            id='cell_id_checklist',
            values=['A'],
        ),

        dcc.Checklist(
            options=[
                {'label': 1, 'value': 1},
                {'label': 2, 'value': 2},
                {'label': 3, 'value': 3},
                {'label': 4, 'value': 4},
                {'label': 5, 'value': 5},
                {'label': 6, 'value': 6},
            ],
            id='replicate_checklist',
            values=range(1, 7),
        ),


        html.Div([
            dcc.Checklist(
                options=[
                    {'label': 0, 'value': 0},
                    {'label': 96, 'value': 96}
                ],
                id='baseline_time_checklist',
                values=[0, 96],
            ),

            dcc.Checklist(
                options=[
                    {'label': 0.5, 'value': 0.5},
                    {'label': 1, 'value': 1},
                    {'label': 2, 'value': 2},
                    {'label': 3, 'value': 3},
                    {'label': 4, 'value': 4},
                    {'label': 8, 'value': 8},
                    {'label': 12, 'value': 12},
                    {'label': 24, 'value': 24},
                    {'label': 48, 'value': 48},
                    {'label': 72, 'value': 72},
                    {'label': 96, 'value': 96}
                ],
                id='time_checklist',
                values=[0.5, 1, 2, 3, 4, 8, 12, 24, 48, 72, 96]
            )
        ], style={'display': 'block-inline'}),

        dcc.Checklist(
            options=[{'label': g, 'value': g} for g in genes],
            id='gene_checklist',
            values=genes,
        ),

        dcc.RadioItems(
            options=[{'label': 'Normalized', 'value': True}],
            id='norm',
            value=True,
        ),
        dcc.Graph(id='graph'),
    ],
)

@app.callback(
    Output('graph', 'figure'),
    [Input('treatment_checklist', 'values'),
     Input('cell_id_checklist', 'values'),
     Input('time_checklist', 'values'),
     Input('replicate_checklist', 'values'),
     Input('baseline_time_checklist', 'values'),
     Input('gene_checklist', 'values'),
     Input('norm', 'value')]
)
def plot_graph(treatment, cell_id,
               time, replicate, baseline_time, gene, norm):

    print treatment
    print cell_id
    print time
    print baseline_time
    print gene
    print norm
    if not isinstance(treatment, list):
        treatment = treatment

    for treat in treatment:
        print treat, type(treat)
        if treat == u'Control':
            control = Query(exp, 'Control', cell_id=cell_id, time=time,
                 replicate=replicate, gene=gene, normed=norm).result.loc['Control']
        elif treat == u'TGFb':
            tgfb = Query(exp, 'TGFb', cell_id=cell_id, time=time,
                            replicate=replicate, gene=gene, normed=norm).result.loc['TGFb']
        elif treat == u'Baseline':
            baseline = Query(exp, 'Baseline', cell_id=cell_id, time=baseline_time,
                            replicate=replicate, gene=gene, normed=norm).result.loc['Baseline']

    print control

    return {
        'data': [go.Scatter(
            x=range(12),
            y=range(12),
        )],
    }



# print plot_graph()



if __name__ == '__main__':
    app.run_server(debug=True)

























