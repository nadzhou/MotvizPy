

import numpy as np

def seq_extract(in_file): 
    seqs = []
    
    with open(in_file, "r") as f: 
        for line in f: 
            line = line.rstrip("\n")
            if ">" not in line: 
                seqs.append(line)
    
    return seqs
class Alignment: 
    def __init__(self):
        self.seq = seq
    
    def seq2np(seq): 
        return np.array(seq)    

    
ar = seq_extract("aligned1.fasta")
np_ar = seq2np(ar)

print(np_ar)