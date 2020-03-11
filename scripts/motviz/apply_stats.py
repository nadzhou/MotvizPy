

import numpy as np

def seq_extract(in_file): 
    seqs = []
    
    with open(in_file, "r") as f: 
        seqs = []
        for line in f: 
            line = line.rstrip("\n")
            if ">" not in line: 
                seqs.append(np.array(list(line), dtype=object))    
        seqs = [[i for i in j] for j in seqs]
        return seqs


def vert_seq_parser(np_seq):    
    whole_seq = np.array(np_seq, dtype=object)
    print(whole_seq.ndim)    
    def _vert_seq(array): 

        return array
    r =  np.apply_along_axis(_vert_seq, 0, whole_seq)   
    print(r[0])

    
ar = seq_extract("/home/nadzhou/blastdb/refseq_protein.00/aligned1.fasta")

vert_seq_parser(ar)
