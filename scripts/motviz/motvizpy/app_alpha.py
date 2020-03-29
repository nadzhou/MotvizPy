import base64
import datetime
import io

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table
from stats import seq_extract
from stats import Analysis

import pandas as pd
import numpy as np

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = (
        html.Div([
            html.H1("MotvizPy - A Motif Discovery tool", 
                    style={'textAlign' : 'center'}),
            
            html.H4("A product of Exon", 
                    style={'textAlign' : 'center'}),
            
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
)


def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)

    try:
        if 'fasta' in filename:
            # Assume that the user uploaded a CSV file
            seq = seq_extract(io.StringIO(decoded.decode('utf-8')))
            seq = [[x for x in y] for y in seq]
            print(seq)
            c = Analysis(seq, "1xef")
            c_ent = c.conservation_score(c.seq2np())
            norm_data = c.normalize_data(c_ent)
            norm_data_len = np.arange(len(norm_data))
            
            print(norm_data)
            
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
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
                        {'x': norm_data_len, 'y': norm_data, 'type': 'line', 'name': 'SF'},
                    ],
                    'layout': dict(
                        xaxis={'title': 'Amino acid position'},
                        yaxis={'title': 'Conservation score'},
                        hovermode='Closest'
                    )
                },
            )
        ])
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
