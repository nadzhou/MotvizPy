

import numpy as np

def seq_extract(in_file): 
    seqs = []
    
    with open(in_file, "r") as f: 
 
        for line in f: 
            line = line.rstrip("\n")
            if ">" not in line: 
                seqs.append(list(line))        
        return seqs

def seq2np(seq): 
    return np.array(seq, dtype=object)

def vert_seq_parser(np_seq):     
    def _vert_seq(array): 

        return array
    r =  np.apply_along_axis(_vert_seq, 0, np_seq)   
    print(r)

    
ar = seq_extract("/home/nadzhou/blastdb/refseq_protein.00/aligned1.fasta")

s = seq2np(ar)

vert_seq_parser(s)

