# -*- coding: utf8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from parse import *
import functools32


app = dash.Dash()

app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

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


GRAPH_WIDTH = 600
GRAPH_HEIGHT = 400

markdown1 = """
TGFb is involved in the maintenance of dermal extracellular matrix (ECM). 
Changes associated with the dermal ECM with age contribute towards
the ageing phenotype in skin tissue. This experiment was designed 
to investigate how fibroblast respond to TGFb over time and how the measured
components differ between three classes of cell:

* Neonatal Human Dermal Fibroblasts (HDFn). ID's = ['A', 'B', 'C']
* Senescent Human Dermal Fibroblasts (Sen). IDs = ['D', 'E', 'F']
* Adult Human Dermal Fibroblasts (A). IDs = ['G', 'H', 'I']

The cells used in the senescent group were the same cell type as the 
neonatal group but senescence was induce with x-ray irradiation.
"""

markdown2 = """
Fibroblasts were treated with 5ng/mL TGFb or PBS (control) for 0.5, 1, 2, 3, 
4, 8, 12, 24, 48, 72 or 96h (n=6). Both control and TGFb time courses were 
serum starved for 24h prior to experiment to remove residual TGFb. 
This is in contrast to the Baseline groups which were collected at 0h and 
96h time points and were not starved or treated in any way. 
"""

markdown3 = """
The entire experiment was conducted in triplicate (9 cell lines total). The 
three sub-experiments each consist of a HDFn, Sen and Adult cell line. 

1) A, D, G
2) B, E, H
3) C, F, I

"""

app.layout = html.Div(
    style={'textAlign': 'center'},
    children=[

        ## Markdown Text
        html.Div([html.H1('Analysis of WaferGen Data for Large Study', style={'textAlign': 'center'}),
                  html.H2('Background', style={'textAlign': 'center'}),
                  dcc.Markdown(markdown1),
                  html.H2('Treatment', style={'textAlign': 'center'}),
                  dcc.Markdown(markdown2),
                  html.H2('Sub Experiments', style={'textAlign': 'center'}),
                  dcc.Markdown(markdown3),
                  html.H1('Options', style={'textAlign': 'center'})
                  ],
                 style={
                     'display': 'inline-block',
                     'textAlign': 'justify',
                     'width': 700
                 }),

        html.Div(
            style={'diplay': 'inline'},
            children=[
            html.Div([
                html.H3('Treatments'),
                dcc.Checklist(
                    options=[
                        {'label': 'TGFb', 'value': 'TGFb'},
                        {'label': 'Control', 'value': 'Control'},
                        {'label': 'Baseline', 'value': 'Baseline'},
                    ],
                    id='treatment_checklist',
                    values=['TGFb', 'Control', 'Baseline'],
                    labelStyle={'display': 'inline-block'}

                )],
                ),

            html.Div([
                html.H3('Replicates'),
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
                    labelStyle={'display': 'inline-block'}

                )],
            ),

            html.Div([
                html.H3('Baseline Time Points'),
                dcc.Checklist(
                    options=[
                        {'label': 0, 'value': 0},
                        {'label': 96, 'value': 96}
                    ],
                    id='baseline_time_checklist',
                    values=[0, 96],
                    labelStyle={'display': 'inline-block'}
                )
                ],
            ),

            html.Div([
                html.H3('Time Points'),
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
                    values=[0.5, 1, 2, 3, 4, 8, 12, 24, 48, 72, 96],
                    labelStyle={'display': 'inline-block'}
                )],
            ),

            html.Div([
                html.H3('Raw Data (Ct) or Normalized to Geometic Mean'),
                dcc.RadioItems(
                    options=[{'label': 'Normalized', 'value': 'Norm2ref'},
                             {'label': 'Raw', 'value': 'Ct'}],
                    id='norm',
                    value='Norm2ref',
                    labelStyle={'display': 'inline-block'}
                )],
            ),

            html.Div([
                html.H3('Set Y axis limit to maximum for gene'),
                dcc.RadioItems(
                    options=[{'label': 'True', 'value': True},
                             {'label': 'False', 'value': False}],
                    id='ylim',
                    value=True,
                    labelStyle={'display': 'inline-block'}
                )],
            ),

            html.Div([
                html.H3('Secondary Normalization'),
                dcc.RadioItems(
                    options=[{'label': 'None', 'value': None},
                             {'label': 'Fold Change', 'value': 'fold_change'}],
                    id='secondary_norm',
                    value=None,
                    labelStyle={'display': 'inline-block'}
                )],
            ),


                html.Div([
                html.H2('Genes to Plot'),
                dcc.Dropdown(
                    options=[{'label': g, 'value': g} for g in genes],
                    id='gene_dropdown',
                    value=genes[0])
            ])
            ],
        ),

        html.Div([
            dcc.Graph(id='A:Neonatal1'),
            dcc.Graph(id='B:Neonatal2'),
            dcc.Graph(id='C:Neonatal3')
        ], style={'float': 'left',
                  'width': GRAPH_WIDTH}),

        html.Div([
            dcc.Graph(id='D:Senescent1'),
            dcc.Graph(id='E:Senescent2'),
            dcc.Graph(id='F:Senescent3')
        ], style={'float': 'left',
                  'width': GRAPH_WIDTH}),

        html.Div([
            dcc.Graph(id='G:Adult1'),
            dcc.Graph(id='H:Adult2'),
            dcc.Graph(id='I:Adult3')
        ], style={'float': 'left',
                  'width': GRAPH_WIDTH})
    ])


def get_highest(gene_name, norm):
    """
    Get highest reading for gene_name accross all
    cell lines, treatments and time points
    :param gene_name:
    :return:
    """
    if norm not in ['Ct', 'Norm2ref']:
        raise ValueError

    if gene_name not in genes:
        raise ValueError('"{}" not in "{}"'.format(gene_name, genes))

    d = [i for i in exp.treatment_data[norm].groupby(level=[3])]
    for i in d:
        if i[0] == gene_name:
            return i[1].max().max()

# print get_highest('SMAD7', 'Norm2ref')
def graphs(graph_id):
    @app.callback(
        Output(graph_id, 'figure'),
        [Input('treatment_checklist', 'values'),
         Input('time_checklist', 'values'),
         Input('replicate_checklist', 'values'),
         Input('baseline_time_checklist', 'values'),
         Input('gene_dropdown', 'value'),
         Input('norm', 'value'),
         Input('ylim', 'value'),
         Input('secondary_norm', 'value')]
    )
    def plot_graph(treatment, time, replicate, baseline_time, gene, norm, ylim,
                   secondary_norm):
        print 'treatments:', treatment
        print 'time:', time
        print 'baseline_time:', baseline_time
        print 'gene:', gene, type(gene)
        print 'norm:', norm
        print 'replicate:', replicate
        print 'ylim', ylim
        if not isinstance(treatment, list):
            treatment = [treatment]

        if not isinstance(gene, list):
            gene = [gene]

        cell_id = graph_id.split(':')[0]

        treatment_data = exp.treatment_data.loc[cell_id][norm]
        baseline_data = exp.baseline_data.loc[cell_id][norm]
        control_data = ''
        tgf_data = ''
        base_data = ''

        if secondary_norm is None:
            data = []
            for treat in treatment:
                for rep in replicate:
                    for g in gene:
                        if treat == u'Control':
                            control_data = treatment_data.loc[treat, rep, g].values
                            data.append(go.Scatter(x=exp.time, y=control_data,
                                        name='{}_{}_{}'.format(treat, rep, g)))

                        elif treat == u'TGFb':
                            tgf_data = treatment_data.loc[treat, rep, g].values
                            data.append(go.Scatter(x=exp.time, y=tgf_data,
                                        name='{}_{}_{}'.format(treat, rep, g)))

                        elif treat == u'Baseline':
                            base_data = baseline_data.loc[treat, rep, g].values
                            data.append(go.Scatter(x=[0, 96], y=base_data,
                                        name='{}_{}_{}'.format(treat, rep, g)))
        elif secondary_norm is 'fold_change':
            pass
        
        highest = max([get_highest(i, norm) for i in gene])
        if ylim:
            print 'ylim is', ylim
            layout = go.Layout(title=graph_id,
                               xaxis=dict(title='Time(h)'),
                               yaxis=dict(title='AU',
                                          range=[0, highest]))
        else:
            print 'ylim not'
            layout = go.Layout(title=graph_id,
                               xaxis=dict(title='Time(h)'),
                               yaxis=dict(title='AU'))

        return go.Figure(data=data, layout=layout)

graph_ids = ['A:Neonatal1', 'B:Neonatal2', 'C:Neonatal3',
             'D:Senescent1', 'E:Senescent2', 'F:Senescent3',
             'G:Adult1', 'H:Adult2', 'I:Adult3']

for g in graph_ids:
    graphs(g)




# print plot_graph()



if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', use_reloader=True, port=12345)




























