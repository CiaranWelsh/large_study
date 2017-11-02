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

























