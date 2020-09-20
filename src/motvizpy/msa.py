#!/usr/bin/env python3.7

from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess

def msa(in_file, out_file): 
    """Take the sequences from the PSI-BLAST result and
    do a Clustal Omega alignment on all the sequences. 
    
    Args: 
        in_file [str]: Path to the input file. 
        
    Returns: 
        out_file [file]: Aligned sequences that is ready for analysis. 
    
    """
        
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, \
                        outfile=out_file, verbose=True, force=True)

    cmd_str = str(clustalomega_cline).split(" ")
    cmd_str[0] = "./c"
    
    return subprocess.run(cmd_str, check=True)



def main(): 

    in_file = "/home/nadzhou/Desktop/divided/fast00.fasta"
    out_file = "/home/nadzhou/Desktop/out.fasta"

    msa(in_file, out_file)





if __name__  == '__main__': 
    main()



