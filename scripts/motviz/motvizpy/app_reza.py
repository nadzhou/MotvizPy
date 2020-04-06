import plotly.graph_objects as go

import numpy as np
import pandas as pd
from pathlib import Path
from plotly.subplots import make_subplots

def view_graph(file, file2): 
    """Visualiza the conservation score dataa from directory
    Graph the data as an interative scatter plot. 
    
    Args: 
        file [str]: Path address for the csv file 
    
    """
    df = pd.read_csv(file)
    df2 = pd.read_csv(file2)
    
    norm_data_len = np.array(df['Nucleotide position'])
    norm_data = np.array(df['Variation score'])

    aa = np.array(df2['Physiochemical properties'])
    freq = np.array(df2["Percent frequency"])
    
    data1 = dict(
        x=norm_data_len, 
        y=norm_data,
        name="Entropy"
    )
    
    data2 = dict(
        labels=aa, 
        values=freq, 
        name="Structure physiochemical properties"
    )
    
    data = [data1, data2]
    
    fig = make_subplots(rows=1, cols=2, specs=[[{'type':'scatter'}, {'type':'pie'}]])
    fig.add_trace(go.Pie(labels=aa, values=freq, scalegroup=1), 1, 2)
    fig.add_trace(go.Scatter(x=norm_data_len, y=norm_data, name="Entropy"), 1, 1)
    fig.update_layout(title_text="Amino Acid positions",
                  title_font_size=30)
    fig.show()  
    

def lay(): 
    """Layout design for the Plotly object. 
    Labels the x and y axes and gives hover animations 
    
    Returns: 
        layout [plotly go object]: Layout object
    
    """
    
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
    path = Path( "/home/nadzhou/Desktop/results.csv")
    path2 = Path("/home/nadzhou/Desktop/plot.csv")
    if path: 
        view_graph(path, path2)
    
    
if __name__ == "__main__": 
    main()