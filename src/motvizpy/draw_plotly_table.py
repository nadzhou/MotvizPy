import pandas as pd
import plotly.graph_objects as go



def draw_table(file_path): 
    df = pd.read_csv(file_path)

    fig = go.Figure(data=[go.Table(
    header=dict(values=list(df.columns),
                fill_color='paleturquoise',
                align='left'),
    cells=dict(values=[df.pdb_id, df.expmethod, df.resolution, df.nr_residues, df.nr_atoms, df.status, df.residues],
            line_color='darkslategray',            
            fill=dict(color=['purple', 'white']),
            align=['left', 'center'],            
            font_size=12,
            height=30))
            
])

    fig.show()








def main(): 
    file_path = "/home/nadzhou/Desktop/pdb_results.csv"

    draw_table(file_path)


if __name__ == '__main__': 
    main()

