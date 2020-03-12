
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
        print(array)
        return array
    return np.apply_along_axis(_per_column, 0, np_seq)
    

    
    
seq = seq_extract("/home/nadzhou/Desktop/aligned1.fasta")

seq = [list(x) for x in seq]
np_seq = np.asarray(seq, dtype='S1')

v = vert_seq(np_seq)

# np_seq = np.transpose(np_seq)

# for i in np_seq: 
#     print(i)

# Started something complete different from here. 
# Truing out the vert_seq function
# seems to be working fine with hard-coded lists
# but has a stick up its ass when I try to do 
# list comprehension. Space. 

# s1 = ["happy", "birthday"]
# s1 = [[x for x in y] for y in s1]

# twod_list = [["1","2","3"],["4","5","6"],["7","8","9"]]
# s2 = vert_seq(s1)
# print(s2.ndim)

# two = vert_seq(twod_list)
# print(two)

# chars_in_y_axis = []
#np.apply_along_axis(lambda x: v_seq.append(x), 0, np_seq)

