from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import seaborn as sns 
import matplotlib.pyplot as plt
import pandas as pd

basic_aa = ['R' , 'H' , 'K' ]
acidic_aa = ['D' , 'E']
hydroxylic_aa = ['S' , 'T' ]
amidic_aa = ['N' , 'Q' ]
aliphatic_aa = ['A' , 'G' , 'L' , 'I' , 'P' , 'V' ]
aromatic_aa = ['F' , 'W' , 'Y' ]

def count_prop(aa_count, prop): 
    """Check the amino acid dict for physiochemical properties 
    
    Args: 
        aa_count [dict]: Dictionary of amino acids with counts for each 
        prop [list: Amino acids from acidic, basic etc list
        
    Returns: 
        prop_percent [list]: Total percentage of the group of amino acids 
    """
    
    prop_percent = [v for i, (k, v) 
            in enumerate(aa_count.items()) if k in prop]
    
    return(sum(prop_percent) * 100)

def seq_properties(file_path): 
    """Apply protein analysis on a fasta file to get analyzed amino acid profile 
    
    Args: 
        file_path [str]: File directory for the fasta file
        
    Returns: 
        total_percent_dict [dict]: Amino acid with counts dict
    """
    
    record = SeqIO.read(file_path, 'fasta')
    analyzed_seq = ProteinAnalysis(str(record.seq))
    
    c = analyzed_seq.get_amino_acids_percent()
    
    acidic_percent = count_prop(c, acidic_aa)
    basic_percent = count_prop(c, basic_aa)
    hydroxylic_percent = count_prop(c, hydroxylic_aa)
    amidic_percent = count_prop(c, amidic_aa)
    aliphatic_percent = count_prop(c, aliphatic_aa)
    aromatic_percent = count_prop(c, aromatic_aa)
    
    total_percent_dict = {
        "Acidic" : acidic_percent, 
        "Basic" : basic_percent, 
        "Hydroxilic" : hydroxylic_percent, 
        "Amidic" : amidic_percent, 
        "Aliphatic" : aliphatic_percent, 
        "Aromatic" : aromatic_percent
        }
    return total_percent_dict
    
def csv_writer(data, file): 
    """Write the data in dict form to a csv file for plotting 
    
    Args: 
        data [dict]: Amino acid char with counts for each
        file [str]: Path for where the amini acid data should be written 
    
    """
    
    df = pd.DataFrame.from_dict(data, orient="index")
    df.index.name = "Physiochemical properties"
    df.columns = ["Percent frequency"]
    df.to_csv(file)
     
def main(): 
    file = '/home/nadzhou/Desktop/1yu5.fasta'
    csv_file = "/home/nadzhou/Desktop/plot.csv"

    c = seq_properties(file)  
    print(c)
    csv_writer(c, csv_file)  

if __name__ == '__main__': 
    main()
    
    
    