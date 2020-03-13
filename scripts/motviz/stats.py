
import numpy as np
from Bio import SeqIO
import math
import matplotlib.pyplot as plt
import seaborn as sns

def seq_extract(in_file): 
    """Extrct the sequence from the fasta file
    Args: 
        input file [str]
    REturns: 
        Sequence [1d list]
    """
    record = SeqIO.parse(in_file, "fasta")
    seqs = []
    for i in record: 
        seqs.append(str(i.seq))
    return seqs
class Analysis: 
    """Class to analyze the multiple sequence alignment. 
    Implements the Shannon Entropy, chi-squared and 
    per-column variation analysis algorithms. 
    """        
    def __init__(self, seq): 
        """
        Initialize class Analysis. 
        
        Args: 
            Sequence list [list of lists]
            
        Returns: 
            Conservation score 
            Chi-squared table 
            etc etc
        """
        def seq2np(seq): 
            return np.asarray(seq, dtype=object)
        
        np_seq = seq2np(seq)
        self.seq = seq
        self.np_seq = np_seq   
    

    def conservation_score(self): 
        """
        Calculate the Shannon Entropy vertically
        for each position in the amino acid msa sequence.
        
        Args: 
            Numpy sequence [nd array]
        Returns: 
            Vertical-orientation conservation score [nd array]    
        """

        def _shannon(array): 
            array = str(array)
            entropy = 0 
            for x in array: 
                char_count = np.char.count(array, sub=x)

                if x == "-": 
                    p_x = float(char_count/len(array)*2) 
                    
                else: p_x = float(char_count/len(array))
                if p_x > 0: 
                    entropy += (- p_x*math.log(p_x, 2)/20)
            return entropy
        return np.apply_along_axis(_shannon, 0, self.np_seq)      

    def normalize_data(self, ent_list): 
        """Takes the entropy array and normalizes the data. 
        
        Args: 
            Shannon Entropy array [numpy ndarray]
        Returns: 
            Normalized numpy array [numpy ndarray]
        """
        return (2.*(ent_list - \
            np.min(ent_list))/np.ptp(ent_list)-1)       
        

        
        
seq = seq_extract("/home/nadzhou/Desktop/aligned1.fasta")
seq = [[x for x in y] for y in seq]

c = Analysis(seq)

c_ent = c.conservation_score()
a = c.normalize_data(c_ent)
a_len = [x for x in range(len(a))]
sns.lineplot(x=a_len, y=a)
plt.show()