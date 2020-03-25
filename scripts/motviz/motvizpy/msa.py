#!/usr/bin/env python3.7

from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess

def msa(in_file): 
    """Take the sequences from the PSI-BLAST result and
    do a Clustal Omega alignment on all the sequences. 
    
    Args: 
        in_file [str]: Path to the input file. 
        
    Returns: 
        out_file [file]: Aligned sequences that is ready for analysis. 
    
    """
        
    out_file = str("/home/nadzhou/Desktop/clustal.fasta")

    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, \
                        outfile=out_file, verbose=True)

    cmd_str = str(clustalomega_cline).split(" ")
    cmd_str[0] = "./c"
    print(cmd_str)

    r = subprocess.Popen(cmd_str)
    if (r.communicate()): 
        return (f"\n Clustal Omega done\nFile written at {out_file}\n")



