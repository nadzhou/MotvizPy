from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess

def msa(in_file): 
    """Generate an MSA of the PSI-BLAST results
    Args 
        input psi_blast file
    Returns 
        Fasta file"""

        
    out_file = "aligned1.fasta"

    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, \
                        outfile=out_file, verbose=True)

    cmd_str = str(clustalomega_cline).split(" ")
    cmd_str[0] = "./c"
    print(cmd_str)

    r = subprocess.Popen(cmd_str)
    if (r.communicate()): 
        print("\n Clustal Omega done")
        print(f"File written at: {out_file} ")


msa("navid.fasta")