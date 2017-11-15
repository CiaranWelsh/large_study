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

design_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'GSS2375_WB_NewDur_Grant')
design_file = os.path.join(design_file, 'new_design.csv')
if not os.path.isfile(design_file):
    os.chdir('..')
    design_file = os.path.abspath(design_file)

if not os.path.isfile(design_file):
    raise ValueError('"{}" not valid file'.format(design_file))

# design_file = os.path
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
                html.H3('Set Y axis limit to'),
                dcc.RadioItems(
                    options=[{'label': 'Max', 'value': 'to_max'},
                             {'label': 'False', 'value':  False}],
                    id='ylim',
                    value=False,
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
                    value=genes[0])]
            )],
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
        print 'secondary_norm', secondary_norm
        time = sorted(time)
        if not isinstance(treatment, list):
            treatment = [treatment]

        if not isinstance(gene, list):
            gene = [gene]

        cell_id = graph_id.split(':')[0]

        treatment_data = exp.treatment_data.loc[cell_id][norm][time]
        baseline_data = exp.baseline_data.loc[cell_id][norm][baseline_time]
        control_data = {}
        tgf_data = {}
        base_data = {}

        highest = {}
        lowest = {}
        d = [i for i in exp.treatment_data[norm][time].groupby(level=[3])]
        for i in d:
            for g in gene:
                if i[0] == g:
                    highest[g] = i[1].max().max() + 0.1 * i[1].max().max()
                    lowest[g] = i[1].min().min() + 0.1 * i[1].min().min()

        if highest == {}:
            raise ValueError('highest is {}')

        if lowest is {}:
            raise ValueError

        colour = 'Pastle2'
        data = []
        if secondary_norm is None:
            for treat in treatment:
                for rep in replicate:
                    for g in gene:
                        if treat == u'Control':
                            control_data[g] = treatment_data.loc[treat, rep, g].values
                            data.append(go.Scatter(x=time,
                                                   y=control_data[g],
                                                   name='{}_{}_{}'.format(treat, rep, g),
                                                   line={'color': colour}
                                                   )
                                        )

                        elif treat == u'TGFb':
                            tgf_data[g] = treatment_data.loc[treat, rep, g].values
                            data.append(go.Scatter(x=time, y=tgf_data[g],
                                                   name='{}_{}_{}'.format(treat, rep, g),
                                                   line={'color': colour}
                                                   )
                                        )

                        elif treat == u'Baseline':
                            base_data[g] = baseline_data.loc[treat, rep, g].values
                            data.append(go.Scatter(x=[0, 96], y=base_data[g],
                                                   name='{}_{}_{}'.format(treat, rep, g),
                                                   line={'color': colour}
                                                   )
                                        )

            if ylim == 'to_max':
                layout = go.Layout(title=graph_id,
                                   xaxis=dict(title='Time(h)'),
                                   yaxis=dict(title='AU',
                                              range=[0, highest[g]]))

            elif ylim == False:
                layout = go.Layout(title=graph_id,
                                   xaxis=dict(title='Time(h)'),
                                   yaxis=dict(title='AU'))\

            return go.Figure(data=data, layout=layout)

        elif secondary_norm == 'fold_change':
            for treat in treatment:
                for rep in replicate:
                    for g in gene:
                        control_data = treatment_data.loc['Control', rep, g].values
                        tgf_data = treatment_data.loc['TGFb', rep, g].values
                        data.append(go.Scatter(x=time, y=tgf_data/control_data,
                                               name='{}_{}_{}'.format(treat, rep, g)))
            print 'ylim --> ', ylim 
            if ylim == 'to_max':
                layout = go.Layout(title=graph_id,
                                   xaxis=dict(title='Time(h)'),
                                   yaxis=dict(title='TGFb / Control (Fold Change)'))

            elif ylim == False:
                layout = go.Layout(title=graph_id,
                                   xaxis=dict(title='Time(h)'),
                                   yaxis=dict(title='TGFb / Control (Fold Change)'))

            return go.Figure(data=data, layout=layout)

graph_ids = ['A:Neonatal1', 'B:Neonatal2', 'C:Neonatal3',
             'D:Senescent1', 'E:Senescent2', 'F:Senescent3',
             'G:Adult1', 'H:Adult2', 'I:Adult3']

for g in graph_ids:
    graphs(g)




# print plot_graph()



if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', use_reloader=True, port=12345)




























