import os
import site
site.addsitedir(os.path.dirname(os.path.dirname(__file__)))
from parse import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import Imputer
import logging
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from flask import Flask
import seaborn
import pandas

logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger()


server = Flask(__name__)

app = dash.Dash(server=server)

app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
#
# from app import app


design_file = os.path.join(
    os.path.dirname(
            os.path.dirname(__file__)), 'GSS2375_WB_NewDur_Grant')

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

factors = sorted(['cell_id', 'treatment', 'replicate', 'sub_experiment',
            'cell_line', 'Treatment Start Date', 'Filename', 'Sample', 'Assay',
            'time_point'])


markdown1 = """
PCA is a dimensionality reduction technique. 
The idea is to reduce high dimensional data sets to 2 or 3
dimensions such that they can be visualized on 2 or 3 
dimensional graphs. The present data set has the following 
variables:

* Cell ID: A-I, n=9
* Cell Line: HDFn, Sen or Adult, n=3
* Treatment: TGFb or Control (Baseline not included here)
* Replicate: 1 - 6 
* Sub-Experiment: n=3
* Treatment start date: Date samples were collected
* Time points: n=11
* Normalized: Boolean
    
Here we lay the data out with time in the columns and all other variables 
in the rows. We then split the data at which ever variable is chosen 
before doing a PCA on the split data and reducing the 11 time points 
to 2D data. Generally the two components explains ~97% of the variance. 
        
PCA cannot handle missing values. Our strategy for handling missing data 
is to set a threshold for the maximum number of missing values per time course.
Time courses with more than this number of missing values are deleted 
before the missing data is imputed using either the mean or median of 
the time course.
"""

app.layout = html.Div([

    html.Div([
        html.H1('Principle Component Analysis'),
        dcc.Link('Go to Time Course Page', href='/apps/plot_web'),
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


    html.Div([
        dcc.Graph(
            id='pca_graph',
            style={'height': 1000}
        ),
    ]),
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
        )],

        dcc.RadioItems(
            id='query-radio',
            options=[
                {'label': "treatment == 'TGFb'", 'value': "treatment == 'TGFb'"}
            ]
        ),

        # dcc.Input(id='input-box'),
        # dcc.Graph(id='pca'),
        #
        # html.Div(id='query', style={'display': 'none'}),
        # #
        # dcc.Input(id='query-box'),

        # html.Button(id='submit-button', n_clicks=0, children='submit'),
    )],
    style={
        'textAlign': 'center',
        'margin-left': 'auto',
        'margin-right': 'auto',
    }

)


# @app.callback(
#     Output('factor', 'multi'),
#     [Input('multi_dropdown_bool', 'value')]
# )
# def set_multi_dropdown(boolean):
#     return boolean

"""
Build query interface. 

List of common queries and an input box to create your own
"""

exp = Experiment(design_file)


# def do_pca_groupby(thresh, strategy, norm, factor=None,
#                    saturation='50%', lightness='50%'):
# @app.callback(
#     Output('query', 'children'),
#     [Input('query-box', 'value')]
# )
# def create_query(query):
#     return "{}".format(query)

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
    treatment_data = exp.treatment_data[norm].stack()
    baseline_data = exp.baseline_data[norm].stack()
    data = pandas.concat([treatment_data, baseline_data])
    data = pandas.DataFrame(data)
    data.columns = [norm]
    '''
    How many measurements shold I have?

    1296 samples, each with 72 readings --> 93312
    Many have been removed automatically due to poor quality
    '''
    data = data.reset_index()

    data = data.rename(columns={
        'cell_line': 'cell_id',
        # 'time_point': 'time',
    })
    # # print treatment.head()
    design = exp.design
    design = design.reset_index()
    print design.keys()    # print design.head()
    data = data.merge(design, on=['cell_id', 'replicate', 'treatment'])

    data = data.drop([u'Lot.Number', u'WG.Plate', u'factor1', u'factor2',
                      u'Iso.Order', u'Iso.Row', u'Label.Row', u'Label.Order',
                      u'Label.Col', u'Iso.Col', u'Batch', u'RNAYield(ug)'], axis=1)
    #'cell_id', 'treatment', 'replicate', 'sub_experiment',
            # 'cell_line', 'Treatment Start Date', 'Filename', 'Sample', 'Assay',
            # 'time_point'
    print data.keys()
    data = data.pivot_table(columns=['Sample', 'cell_id', 'replicate',
                                     'treatment', 'time_point', 'Treatment Start Date',
                                     'sub_experiment', 'Filename', 'cell_line'], index='Assay')
    data = data.drop('time', axis=1)
    if data.shape != (72, 1296):
        raise ValueError

    # pca = PCA(2)
    # pca.fit(data)
    # print pca.explained_variance_
    # data = data.set_index(factors)
    data.dropna(axis=0, how='all', thresh=thresh, inplace=True)
    I = Imputer(axis=1, strategy=strategy)
    data = pandas.DataFrame(I.fit_transform(data), index=data.index, columns=data.columns)
    data = data.transpose()
    pca = PCA(10)
    pca.fit(data)

    # pandas.DataFrame(pca.explained_variance_).plot()
    # plt.show()

    df = pandas.DataFrame(pca.transform(data), index=data.index)
    df = df[[0, 1]]

    return {
        'data': go.Scatter(x=range(10), y=range(10)),
        'layout': go.Layout(title='t')
    }
    #.query("treatment == 'TGFb'")
    # colors = ['hsl({},{}%,{}%)'.format(h, saturation, lightness) for h in numpy.linspace(0, 300, len(groupby_obj))]
    # labels = []
    # for label, df in groupby_obj:
    #     labels.append(label)
    #
    # colours = dict(zip(labels, colors))
    #
    # traces = []
    # for label, df in groupby_obj:
    #     if isinstance(label, tuple):
    #         name = reduce(lambda x, y: '{}_{}'.format(x, y), label)
    #     else:
    #         name = label
    #     df = df.transpose()
    #     pca = PCA(2)
    #     pca.fit(df)
    #     pc = pandas.DataFrame(pca.transform(df), index=df.index)
    #
    #     trace = go.Scatter(
    #         x=pc[0],
    #         y=pc[1],
    #         mode='markers+text+lines',
    #         name=name,
    #         marker=dict(
    #             color=colours[label],
    #             size=15,
    #         ),
    #         text=['{}h'.format(i) for i in pc.index],
    #         textposition='bottom',
    #     )
    #     traces.append(trace)
    #
    # if isinstance(factor, (str, unicode)):
    #     factor = [factor]
    #
    # layout = go.Layout(
    #     title='PCA Data Split by "{}"'.format(reduce(lambda x, y: "{}, {}".format(x, y), factor)),
    #     titlefont=dict(
    #         size=35
    #     ),
    #     hovermode='closest',
    #     legend=dict(
    #         x=-0.15,
    #         y=1,
    #         bgcolor='#E2E2E2'),
    #     hoverlabel=dict(namelength=-1),
    #     xaxis=dict(
    #         title='PC1',
    #         titlefont=dict(
    #             size=20
    #         )
    #
    #         ),
    #     yaxis=dict(
    #         title='PC2',
    #         titlefont=dict(
    #             size=20
    #         )
    #
    #     )
    # )
    # return {
    #     'data': traces,
    #     'layout': layout
    # }



# do_pca_groupby(5, 'median', 'Ct', )

if __name__ == '__main__':
    app.run_server(debug=True)












