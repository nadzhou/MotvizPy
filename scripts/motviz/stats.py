
import numpy as np
from Bio import SeqIO
import math

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
      

    def __init__(self, seq): 
        """
        Take in a 2d sequence and then do 
        a whole bunch of calculations. 
        
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
        print(self.np_seq)
            
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
                p_x = float(char_count/len(array)) 
                if p_x > 0: 
                    entropy += - p_x*math.log(p_x, 2)
            return entropy
        return np.apply_along_axis(_shannon, 0, self.np_seq)
        

        
        
seq = seq_extract("/home/nadzhou/Desktop/aligned1.fasta")
seq = [[x for x in y] for y in seq]

c = Analysis(seq)

print(c.conservation_score())