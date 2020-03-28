import plotly.graph_objects as go

import numpy as np
import pandas as pd

def view_graph(file): 
    df = pd.read_csv(file)
    
    norm_data_len = np.array(df['Amino acid position'])
    norm_data = np.array(df['Conservation score'])

    print(norm_data)
    data = go.Scatter(
        x=norm_data_len, 
        y=norm_data,
        name="Entropy"
    )
    
    fig = go.Figure(data=data, layout=lay())
    fig.show()  
    

def lay(): 
    layout = go.Layout(
    title=go.layout.Title(
        text='Conservation score per amino acid position',
        xref='paper',
        x=0.5,
        font=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
    ),
    xaxis=go.layout.XAxis(
        title=go.layout.xaxis.Title(
            text='Amino acid position',
            font=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    ),
    yaxis=go.layout.YAxis(
        title=go.layout.yaxis.Title(
            text='Conservation score',
            font=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    ), 
    hovermode="closest",
                     updatemenus=[dict(type="buttons",
                                       buttons=[dict(label="Play",
                                                     method="animate",
                                                     args=[None])])]
)
    return layout


def main(): 
    path = "/home/nadzhou/Desktop/biryani.csv"
    view_graph(path)
    
    
if __name__ == "__main__": 
    main()