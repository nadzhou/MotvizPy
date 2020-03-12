
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

def conservation_score(np_seq): 
    """
    Do a vertical manipulation on the array. 
    Take out the sequence vertically
    
    Args: 
        Numpy sequence [nd array]
    Returns: 
        Vertical sequence from each nd array [nd array]    
    """

    def _shannon(array): 
        array = str(array)
        entropy = 0 
        for x in array: 
            char_count = np.char.count(array, sub=x)
            print(char_count)
            p_x = float(char_count/len(array)) 
            if p_x > 0: 
                entropy += - p_x*math.log(p_x, 2)
        return entropy
    return np.apply_along_axis(_shannon, 0, np_seq)
    

    
    
seq = seq_extract("/home/nadzhou/Desktop/aligned1.fasta")

seq = [list(x) for x in seq]
np_seq = np.asarray(seq, dtype=object)
v = conservation_score(np_seq)

