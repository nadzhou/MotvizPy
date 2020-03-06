

import numpy

def seq_extract(in_file): 
    seqs = []
    
    with open(in_file, "r") as f: 
        for line in f: 
            if ">" not in line: 
                seqs.append(line)
    
    for i in seqs: 
        print(i)
    
    
    
seq_extract("aligned1.fasta")