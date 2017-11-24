#-*- coding: utf-8 -*-
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


design_file = os.path.join(
    os.path.dirname(
            os.path.dirname(__file__)), 'GSS2375_WB_NewDur_Grant')

design_file = os.path.join(design_file, 'new_design.csv')
if not os.path.isfile(design_file):
    os.chdir('..')
    design_file = os.path.abspath(design_file)


if not os.path.isfile(design_file):
    raise ValueError('"{}" not valid file'.format(design_file))

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
    
Here we lay the data out with samples along columns and genes down 
 the rows. The PCA is calculated and the first three PCs are 
displayed in 3D. 
        
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

            html.Div([
                html.H2('Choose data point label'),
                dcc.Checklist(
                    id='text',
                    options=[{'label': i, 'value': i} for i in factors],
                    values=['treatment'],
                    labelStyle={'display': 'inline-block'},
                )
            ]),

            html.Div([
                html.H2('Choose which factor to colour plot by'),
                dcc.Checklist(
                    id='colour_by',
                    options=[{'label': i, 'value': i} for i in factors],
                    values=['treatment'],
                    labelStyle={'display': 'inline-block'},
                )
            ])
        ]),

        html.Div([
            html.H2('Query the Data'),
            dcc.Markdown('Construct a [query](https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.query.html) '
                          'for a pandas data frame to subset the data before plotting'),
            dcc.Input(id='input-1-state', type="text", value='All Data'),
            html.Button(id='submit-button', n_clicks=0, children='Submit'),
            html.Div(id='output-state')
        ]),

    html.Div([
        dcc.Markdown('''## Explained Variance
        * PC1 (x) = 52.27%
        * PC2 (y) = 10.87%
        * PC3 (z) = 7.84%
         '''),
        dcc.Graph(
            id='pca_graph',
            style={'height': 1000}
        ),

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
        )

    ]),
        ]),

    html.Label('Note: click on a legend label to toggle it on or off. Double click a legend label to '
               'isolate that trace'),
    ],
    style={
        'textAlign': 'center',
        'margin-left': 'auto',
        'margin-right': 'auto',
    }

)


exp = Experiment(design_file)


@app.callback(Output('output-state', 'children'),
              [Input('submit-button', 'n_clicks')],
              [State('input-1-state', 'value')])
def update_output(n_clicks, input1):
    return input1

def str_join(df, sep, *cols):
    from functools import reduce
    return reduce(lambda x, y: x.astype(str).str.cat(y.astype(str), sep=sep),
                  [df[col] for col in cols])

def print_full(df):
    default = pandas.get_option('display.max_rows')
    pandas.set_option('display.max_rows', df.shape[0])
    print df
    pandas.set_option('display.max_rows', default)


@app.callback(
    Output('pca_graph', 'figure'),
    [
        Input('nan_thresh', 'value'),
        Input('imputation_strategy', 'value'),
        Input('norm', 'value'),
        # Input('factor', 'value'),
        Input('output-state', 'children'),
        Input('text', 'values'),
        Input('colour_by', 'values'),
        Input('saturation', 'value'),
        Input('lightness', 'value'),
    ]
)
def do_pca_groupby(thresh, strategy, norm,
                   output_state, text, colour_by,
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
            'time': 'time_point'
        })
        # # print treatment.head()
        design = exp.design
        design = design.reset_index()
        # print data.keys()
        # print design.keys()
        # data.to_csv('data.csv')
        data = data.merge(design, on=['cell_id', 'replicate', 'treatment', 'time_point'])

        data = data.drop([u'Lot.Number', u'WG.Plate', u'factor1', u'factor2',
                          u'Iso.Order', u'Iso.Row', u'Label.Row', u'Label.Order',
                          u'Label.Col', u'Iso.Col', u'Batch', u'RNAYield(ug)'], axis=1)
        #'cell_id', 'treatment', 'replicate', 'sub_experiment',
                # 'cell_line', 'Treatment Start Date', 'Filename', 'Sample', 'Assay',
                # 'time_point'
        data = data.pivot_table(columns=['Sample', 'cell_id', 'replicate',
                                         'treatment', 'time_point', 'Treatment Start Date',
                                         'sub_experiment', 'Filename', 'cell_line'], index='Assay')
        # if data.shape != (72, 1296):
        #     raise ValueError('"{}" is not (1296, 72)'.format(data.shape))

        data.dropna(axis=0, how='all', thresh=thresh, inplace=True)
        I = Imputer(axis=1, strategy=strategy)
        data = pandas.DataFrame(I.fit_transform(data), index=data.index, columns=data.columns)
        data = data.transpose()
        pca = PCA(10)

        pca.fit(data)

        explained_var = pandas.DataFrame(pca.explained_variance_)
        print explained_var
        plt.show()
        df = pandas.DataFrame(pca.transform(data))
        df.index = data.index
        df = df[[0, 1, 2]]

        df = df.reset_index()
        df['time_point'].astype(float)
        # df['time_point'].astype(float)
        df = df.sort_values(by=text)


        print 'text is' , text

        if text != 'All Data' or text != '':

            try:
                df = df.query(output_state)
            except SyntaxError:
                print 'Query "{}" caused Syntax error'.format(output_state)

        print_full(df)
        print df.shape
        print str_join(df, '_', *text)

        groupby_obj = df.groupby(by=colour_by)
        import colorlover as cl
        print 'len', len(groupby_obj)
        # try:
        #     paired = cl.scales['12']['qual']['Paired']
        #     colours = cl.interp(paired, len(groupby_obj))
        # except ZeroDivisionError:
        colours = ['hsl({},{}%,{}%)'.format(h, saturation, lightness) for h in numpy.linspace(0, 300, len(groupby_obj))]

        labels = []
        for label, df in groupby_obj:
            labels.append(label)

        colours = dict(zip(labels, colours))

        traces = []
        for label, d in groupby_obj:
            if isinstance(label, tuple):
                name = reduce(lambda x, y: '{}_{}'.format(x, y), label)
            else:
                name = label
            trace = go.Scatter3d(
                x=d[0],
                y=d[1],
                z=d[2],
                mode='markers+text',
                text=str_join(d, '_', *text),
                textposition='bottom',
                marker=dict(
                    color=colours[label]
                ),
                name=name
            )
            traces.append(trace)
        layout = go.Layout(
            title='PCA',
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
                title='PC1 Explained Variance {}%'.format(float(explained_var.iloc[0])),
                titlefont=dict(size=20)
                ),

            yaxis=dict(
                title='PC2 Explained Variance {}%'.format(float(explained_var.iloc[1])),
                titlefont=dict(size=20)
            ),

            # zaxis=dict(
            #     title='PC3 Explained Variance {}%'.format(explained_var.iloc[2]),
            #     titlefont=dict(size=20)
            # )
        )


        # print df
        # import plotly
        # plotly.offline.plot(go.Figure(data=[trace], layout=layout))
        return {
            'data': traces,
            'layout': layout,
        }

# do_pca_groupby(5, 'median', 'Ct', )

if __name__ == '__main__':
    app.run_server(debug=True)












