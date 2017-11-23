# -*- coding: utf8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from parse import *
import functools32
from flask import Flask
import dash_auth
from sklearn.preprocessing import Imputer
from sklearn.decomposition import PCA

VALID_USERNAME_PASSWORD_PAIRS = [
    ['Ciaran', 'Welsh',
     'Daryl', 'Shanley',
     'Stefan', 'Pryzborski',
     'Nicola', 'Fullard',
     'Bob', 'Isfort',
     'Charlie', 'Bascom',
     'Ryan', 'Tassef']
]



server = Flask(__name__)

app = dash.Dash(server=server)

auth = dash_auth.BasicAuth(
    app,
    VALID_USERNAME_PASSWORD_PAIRS
)

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


genes = ['ACTA2', 'ADAMTS1', 'ATP6AP1', 'B2M', 'BGN', 'BHLHE40', 'CAV1',
         'CDKN2A', 'CDKN2B', 'COL1A1', 'COL1A2', 'COL2A1', 'COL4A1', 'COL5A1',
         'CTGF', 'DCN', 'EGR2', 'ELN', 'ENG', 'ETS1', 'FBLN1', 'FBN1', 'FGF2',
         'FN1', 'FOSB', 'GADD45B', 'GUSB', 'HAS2', 'ID1', 'IL1A', 'IL1B', 'IL6',
         'ITGA1', 'ITGA2', 'JUN', 'JUNB', 'LARP6', 'LOXL1', 'LOXL2', 'LTBP2',
         'MMP1', 'MMP13', 'MMP14', 'MMP2', 'NOX4', 'PDGFA', 'PMEPA1', 'PPIA',
         'PPP3CA', 'PSMD14', 'PTEN', 'RARA', 'RARG', 'RHOB', 'SERPINE1',
         'SERPINE2', 'SKI', 'SKIL', 'SMAD3', 'SMAD7', 'SPARC', 'TGFB1',
         'TGFBR1', 'TGFBR2', 'THBS1', 'THBS2', 'TIMP1', 'TIMP3', 'TP53BP1',
         'VCAN', 'VEGFA', 'VIM']

factors = sorted(['cell_id', 'treatment', 'replicate', 'sub_experiment',
            'cell_line', 'Treatment Start Date', 'Filename', 'Sample', 'Assay',
            'time_point'])

GRAPH_WIDTH = 500
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
        html.Div([html.H1('WaferGen Data From Dermal Fibroblasts treated with 5ng/mL TGFb', style={'textAlign': 'center'}),
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
                  'width': GRAPH_WIDTH}),

    html.Div([

        html.Div([
            html.H1('Principle Component Analysis'),
            dcc.Markdown(markdown1)],
            style={'textAlign': 'justify',
                       'width': 800,
                       'margin': 'auto'}),

        html.Div([
            html.Div([
                html.H2('Threshold Slider'),
                dcc.Slider(
                    min=1,
                    max=11,
                    step=1,
                    value=5,
                    marks={i: i for i in range(1, 12)},
                    id='nan_thresh'
                )],
                style={'align': 'center',
                       'margin': 'auto',
                       'width': '60%'}
            ),

            html.Br(),

            html.Div([

                html.H2('Imputation Strategy'),

                dcc.RadioItems(
                    options=[
                        {'label': 'mean', 'value': 'mean'},
                        {'label': 'median', 'value': 'median'},
                    ],
                    value='mean',
                    id='imputation_strategy',
                    labelStyle={'display': 'inline-block'},
                ),

                html.Br(),

                html.H2('Raw Ct Values or After Normalization to geomean(PPIA, B2M)'),

                dcc.RadioItems(
                    options=[
                        {'label': 'Raw Ct Values', 'value': 'Ct'},
                        {'label': 'Normalized delta Ct', 'value': 'Norm2ref'}
                    ],
                    value='Ct',
                    id='norm',
                    labelStyle={'display': 'inline-block'},
                ),

                # dcc.RadioItems(
                #     options=[
                #         {'label': 'Raw Ct Values', 'value': 'Ct'},
                #         {'label': 'Normalized delta Ct', 'value': 'Norm2ref'}
                #     ],
                #     value='Ct',
                #     id='plotly_mode',
                #     labelStyle={'display': 'inline-block'},
                # ),

                html.H2('Pick Factor(s)'),

                html.H3('Multi-dropdown'),
                dcc.RadioItems(
                    options=[
                        {'label': 'True', 'value': True},
                        {'label': 'False', 'value': False},
                    ],
                    labelStyle={'display': 'inline-block'},
                    value=False,
                    id='multi_dropdown_bool',
                ),

                html.Label('Note: In "multi-mode" Order affects the colour scale (i.e. try selecting '
                           'treatments then assay or vice versa)'),

                dcc.Dropdown(
                    options=[{'label': factor, 'value': factor} for factor in factors],
                    value=['time_point'],
                    id='factor',
                ),
            ]),
        ]),

        html.Div([
            dcc.Graph(
                id='pca_graph',
                style={'height': 1000}
            ),
        ]),

        html.Label('Note: click on a legend label to toggle it on or off. Double click a legend label to '
                   'isolate that trace'),

        html.Div([
            html.Br(),
            html.H2('Colour Saturation'),
            dcc.Slider(
                min=0,
                max=100,
                step=1,
                value=50,
                id='saturation'
            ),

            html.Br(),
            html.H2('Colour Brightness'),
            dcc.Slider(
                min=0,
                max=100,
                step=1,
                value=50,
                id='lightness'
            )]
        )],
        style={
            'textAlign': 'center',
            'margin-left': 'auto',
            'margin-right': 'auto',
        }

    )
    ])


@app.callback(
    Output('factor', 'multi'),
    [Input('multi_dropdown_bool', 'value')]
)
def set_multi_dropdown(boolean):
    return boolean



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
            for rep in replicate:
                for g in gene:
                    control_data = treatment_data.loc['Control', rep, g].values
                    tgf_data = treatment_data.loc['TGFb', rep, g].values
                    data.append(go.Scatter(x=time, y=tgf_data/control_data,
                                           name='{}'.format(rep)))
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

@app.callback(
    Output('pca_graph', 'figure'),
    [
        Input('nan_thresh', 'value'),
        Input('imputation_strategy', 'value'),
        Input('norm', 'value'),
        Input('factor', 'value'),
        Input('saturation', 'value'),
        Input('lightness', 'value'),
    ]
)
def do_pca_groupby(thresh, strategy, norm, factor,
                   saturation, lightness):
    """
    Do PCA with the entire dataframe (i.e. one data point per time point

    :return:
    """
    seaborn.set_context(context='talk', font_scale=1.5)
    # print exp.treatment_data.head()
    # print exp.design.head()
    treatment = exp.treatment_data[norm]
    treatment = treatment.reset_index()
    treatment = treatment.rename(columns={
        'cell_line': 'cell_id'
    })
    # print treatment.head()
    design = exp.design
    design = design.reset_index()
    # print design.head()
    data = treatment.merge(design, on=['cell_id', 'replicate', 'treatment'])
    data = data.drop([u'Lot.Number', u'WG.Plate', u'factor1', u'factor2',
                      u'Iso.Order', u'Iso.Row', u'Label.Row', u'Label.Order',
                      u'Label.Col', u'Iso.Col', u'Batch', u'RNAYield(ug)'], axis=1)

    data = data.set_index(factors)
    data.dropna(axis=0, how='all', thresh=thresh, inplace=True)
    I = Imputer(axis=1, strategy=strategy)
    data = pandas.DataFrame(I.fit_transform(data), index=data.index, columns=data.columns)

    groupby_obj = data.groupby(level=factor)

    colors = ['hsl({},{}%,{}%)'.format(h, saturation, lightness) for h in numpy.linspace(0, 300, len(groupby_obj))]
    labels = []
    for label, df in groupby_obj:
        labels.append(label)

    colours = dict(zip(labels, colors))

    traces = []
    for label, df in groupby_obj:
        if isinstance(label, tuple):
            name = reduce(lambda x, y: '{}_{}'.format(x, y), label)
        else:
            name = label
        df = df.transpose()
        pca = PCA(2)
        pca.fit(df)
        pc = pandas.DataFrame(pca.transform(df), index=df.index)

        trace = go.Scatter(
            x=pc[0],
            y=pc[1],
            mode='markers+text+lines',
            name=name,
            marker=dict(
                color=colours[label],
                size=15,
            ),
            text=['{}h'.format(i) for i in pc.index],
            textposition='bottom',
        )
        traces.append(trace)

    if isinstance(factor, (str, unicode)):
        factor = [factor]

    layout = go.Layout(
        title='PCA Data Split by "{}"'.format(reduce(lambda x, y: "{}, {}".format(x, y), factor)),
        titlefont=dict(
            size=35
        ),
        hovermode='closest',
        legend=dict(
            x=-0.15,
            y=1,
            bgcolor='#E2E2E2'),
        hoverlabel=dict(namelength=-1),
        xaxis=dict(
            title='PC1',
            titlefont=dict(
                size=20
            )

            ),
        yaxis=dict(
            title='PC2',
            titlefont=dict(
                size=20
            )

        )
    )
    return {
        'data': traces,
        'layout': layout
    }



# print plot_graph()



if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', use_reloader=True, port=12345)




























