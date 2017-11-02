import dash_core_components as html
import dash_core_components as dcc
import dash
from dash.dependencies import Input, State, Event, Output


## import files as we go for a multipage app


##style sheet. Use the one provided with dash
dcc._js_dist[0]['external_url'] = (
    'https://cdn.plot.ly/plotly-basic-1.31.0.min.js'
)


def create_contents(contents):
    """

    :param contents:
    :return:
    """
    h = []
    for content in contents:
        if isinstance(content, list):
            h.append(create_contents(content))
        else:
            h.append(html.Li(content))
    return html.Ul(h)



# header = html.Div(
#     className='header',
#     children=html.Div(
#         className='container-width',
#         style={'height': '100%'},
#         children=[
#             html.A('https://www.google.co.uk')
#         ]
#     )
# )
#

app = dash.Dash()
server = app.server
app.config.suppress_callback_exceptions = True







