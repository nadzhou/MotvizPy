
import numpy as np
from Bio import SeqIO

def seq_extract(in_file): 
    record = SeqIO.parse(in_file, "fasta")
    seqs = []
    for i in record: 
        seqs.append(i.seq)
        
    return seqs

def seq2np(seqs): 
    return np.array(seq)

def vert_seq(np_seq): 
    def _per_column(array): 
        return array

    return np.apply_along_axis(_per_column, 0, np_seq)
    
seq = seq_extract("/home/nadzhou/Desktop/navid.fasta")


np_seq = np.array([np.array(y, dtype=object) for y in seq], dtype=object)

print(np_seq.ndim)