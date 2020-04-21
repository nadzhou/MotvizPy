import base64
import datetime
import io

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table

import pandas as pd
import numpy as np

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

data = pd.read_csv("/home/nadzhou/Desktop/biryani.csv")
amino_acid_pos = np.array(data["Amino acid position"])
cons_score = np.array(data["Conservation score"])


app.layout = (
    html.Div([
        html.Div(children=[
            html.H1(children='Conservation score per amino acid position', 
                    style={
                        'textAlign': 'center',
                        'font-family' : 'Century Gotic'}
                    ),

            html.Div(children='''
                Interactive representation of conservation score.
            ''', 
            style={
                        'textAlign': 'center',
                        'font-family' : 'Century Gothic'
                    }     
            ),

            dcc.Graph(
                id='example-graph',
                figure={
                    'data': [
                        {'x': amino_acid_pos, 'y': cons_score, 'type': 'line', 'name': 'SF'},
                    ],
                    'layout': dict(
                        xaxis={'title': 'Amino acid position'},
                        yaxis={'title': 'Conservation score'},
                        hovermode='Closest'
                    )
                },
            )
        ]),
        html.Div([
            html.Div([
                dcc.Upload(
                    id='upload-data',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select Files')
                    ]),
                    style={
                        'width': '100%',
                        'height': '60px',
                        'lineHeight': '60px',
                        'borderWidth': '1px',
                        'borderStyle': 'dashed',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                        'margin': '10px'
                    },
                    # Allow multiple files to be uploaded
                    multiple=True
                ),
                html.Div(id='output-data-upload'),
                
            ])
        ])
    ])
)


def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in df.columns]
        ),

        html.Hr(),  # horizontal line

        # For debugging, display the raw contents provided by the web browser
        html.Div('Raw Content'),
        html.Pre(contents[0:200] + '...', style={
            'whiteSpace': 'pre-wrap',
            'wordBreak': 'break-all'
        })
    ])



@app.callback(Output('output-data-upload', 'children'),
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified')])
def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children



    

if __name__ == '__main__':
    app.run_server(debug=True)

# def lay(norm_data, norm_data_len): 
#     html.Div(style={'backgroundColor': colors['background']}, children=[
#     html.H1(
#         children='Conservation Scores',
#         style={
#             'textAlign': 'center',
#             'color': colors['text']
#         }
#     ),

#     html.Div(children='Per base comparison of conservation for amino acid', style={
#         'textAlign': 'center',
#         'color': colors['text']
#     }),

#     dcc.Graph(
#         id='example-graph-2',
#         figure={
#             'data': [
#                 {'x': nor_data_len, 'y': norm_data, 'type': 'bar', 'name': 'SF'},
#             ],
#             'layout': {
#                 'plot_bgcolor': colors['background'],
#                 'paper_bgcolor': colors['background'],
#                 'font': {
#                     'color': colors['text']
#                 }
#             }
#         }