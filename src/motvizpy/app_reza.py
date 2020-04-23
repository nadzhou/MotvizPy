import plotly.graph_objects as go

import numpy as np
import pandas as pd
from pathlib import Path
from plotly.subplots import make_subplots
import plotly.express as px



def view_graph(file): 
    """Visualiza the conservation score dataa from directory
    Graph the data as an interative scatter plot. 
    
    Args: 
        file [str]: Path address for the csv file 
    
    """
    df = pd.read_csv(file)
    
    norm_data_len = np.array(df['Nucleotide position'])
    norm_data = np.array(df['Variation score'])


    fig = go.Figure()

    fig.add_trace(
            go.Scatter(
                x=norm_data_len, 
                y=norm_data,  
            )
    )


    fig.update_layout(lay())


    fig.show()  
    

def draw_pie_chart(file2): 
    df2 = pd.read_csv(file2)

    aa = np.array(df2['Physiochemical properties'])
    freq = np.array(df2["Percent frequency"])

    fig = go.Figure()

    fig.add_trace(
        go.Pie(
            labels=aa, 
            values=freq, 
            hole=.3,
            pull=[0, 0, 0.1, 0],
            name="Structure physiochemical properties", 
        )
    )

    fig.update_layout(uniformtext_minsize=14, uniformtext_mode='hide')


    fig.show()

def draw_table(file_path):
    df = pd.read_csv(file_path)

    fig = go.Figure(data=[go.Table(
    header=dict(values=list(df.columns),
            fill_color='purple',
            font=dict(color='white', size=17),
            align='center'),
    cells=dict(values=[df.pdb_id, df.expmethod, df.resolution, 
                    df.nr_residues, df.nr_atoms, df.status, df.residues],
            line_color='darkslategray',            
            fill=dict(color=['plum', 'white']),
            align='center',            
            font_size=13,
            height=30))
            
])

    fig.show()






def layout2(): 
    layout = go.Layout(
        title=go.layout.Title(
            xref='paper',
            x=0.5,
            font=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                )
        ),
        xaxis = dict(
            title= dict(
                text='Amino acid position',
                font=dict(
                        family='Courier New, monospace',
                        size=18,
                        color='#7f7f7f'
                    )
        )   
    ), 
        yaxis2 = dict(
            title=dict(
                text="Conservation score", 
                font=dict(
                        family='Courier New, monospace',
                        size=18,
                        color='#7f7f7f'
                    )

        )

    )
    

    )
    return layout




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
    file_path = "/home/nadzhou/Desktop/pdb_results.csv"

    if path: 
        view_graph(path)

    draw_pie_chart(path2)
    draw_table(file_path)
    



if __name__ == "__main__": 
    main()