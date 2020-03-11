
import numpy as np
from Bio import SeqIO

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
        seqs.append(i.seq)
        
    return seqs

def vert_seq(np_seq): 
    """
    Do a vertical manipulation on the array. 
    Take out the sequence vertically
    
    Args: 
        Numpy sequence [nd array]
    Returns: 
        Vertical sequence from each nd array [nd array]    
    """
    def _per_column(array): 
        return array

    return np.apply_along_axis(_per_column, 0, np_seq)
    
seq = seq_extract("/home/nadzhou/Desktop/navid.fasta")

# I was trying to do direct conversion 
# but doesn't seem to happen 
np_seq = np.array([np.array(y, dtype=object) for y in seq], dtype=object)

# ndim only gives me 1
# but when I go inside, 
# gives me the internal dim as well 

# possibly numpy problem? 

#print(np_seq)
print(np_seq[0].shape)

ver = vert_seq(np_seq)
print(ver)