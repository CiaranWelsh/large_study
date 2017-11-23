import os
import site
site.addsitedir(os.path.dirname(os.path.dirname(__file__)))
from parse import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import Imputer
import logging
import dash
from dash.dependencies import Input, Output
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


@app.callback(
    Output('factor', 'multi'),
    [Input('multi_dropdown_bool', 'value')]
)
def set_multi_dropdown(boolean):
    return boolean



exp = Experiment(design_file)


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





if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=12345)












