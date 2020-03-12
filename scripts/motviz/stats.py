
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
        seqs.append(str(i.seq))
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
    

    
    
seq = seq_extract("/home/nadzhou/Desktop/aligned1.fasta")

seq = [list(x) for x in seq]
np_seq = np.asarray(seq, dtype='S1')

v = vert_seq(np_seq)
