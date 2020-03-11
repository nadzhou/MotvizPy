

import numpy as np

def seq_extract(in_file): 
    seqs = []
    
    with open(in_file, "r") as f: 
        seqs = []
        for line in f: 
            line = line.rstrip("\n")
            line_list = list(line)
            if ">" not in line: 
                seqs.append(line_list)    
        return seqs
    
def seq2np(seqs): 
    return np.array(seq)

def vert_seq(np_seq): 
    def _per_column(array): 
        # print(f"array {array}")
        trans = array
        # print(f"trans {trans}")
        return trans
    new_ar =  np.apply_along_axis(_per_column, 0, np_seq)
    
    return new_ar
    
seq = seq_extract("/home/nadzhou/Desktop/navid.fasta")

np_seq = np.array(seq, dtype=object)

vert_seq1 = vert_seq(np_seq)

print(vert_seq1)