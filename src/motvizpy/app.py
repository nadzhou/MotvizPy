# import dash
# import dash_core_components as dcc
# import dash_html_components as html

# import numpy as np
# import pandas as pd

# def view_graph(file): 
#     # df = pd.read_csv(file)
    
#     # norm_data_len = np.array([df['Amino acid position']])
#     # norm_data = np.array([df['Conservation value']])

#     # data = go.Scatter(
#     #     x=norm_data_len, 
#     #     y=norm_data,
#     #     name="Entropy"
#     # )
#     # fig = go.Figure(data=data, layout=lay())
#     # fig.show()  
    
#     app = dash.Dash(__name__, external_spreadsheets=file)
    
#     colors = {
#     'background': '#111111',
#     'text': '#7FDBFF'
#     }
    
#     app.layout = 


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
#     )
# ])

# # def lay(): 
# #     layout = go.Layout(
# #     title=go.layout.Title(
# #         text='Conservation score per amino acid position',
# #         xref='paper',
# #         x=0.5,
# #         font=dict(
# #                 family='Courier New, monospace',
# #                 size=18,
# #                 color='#7f7f7f'
# #             )
# #     ),
# #     xaxis=go.layout.XAxis(
# #         title=go.layout.xaxis.Title(
# #             text='Amino acid position',
# #             font=dict(
# #                 family='Courier New, monospace',
# #                 size=18,
# #                 color='#7f7f7f'
# #             )
# #         )
# #     ),
# #     yaxis=go.layout.YAxis(
# #         title=go.layout.yaxis.Title(
# #             text='Conservation score',
# #             font=dict(
# #                 family='Courier New, monospace',
# #                 size=18,
# #                 color='#7f7f7f'
# #             )
# #         )
# #     ), 
# #     hovermode="closest",
# #                      updatemenus=[dict(type="buttons",
# #                                        buttons=[dict(label="Play",
# #                                                      method="animate",
# #                                                      args=[None])])]
# # )
# #     return layout

import dash
import dash_core_components as dcc
import dash_html_components as html

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),

    dcc.Graph(
        id='example-graph',
        figure={
            'data': [
                {'x': [1, 2, 3], 'y': [4, 1, 2], 'type': 'bar', 'name': 'SF'},
                {'x': [1, 2, 3], 'y': [2, 4, 5], 'type': 'bar', 'name': u'Montréal'},
            ],
            'layout': {
                'title': 'Dash Data Visualization'
            }
        }
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)