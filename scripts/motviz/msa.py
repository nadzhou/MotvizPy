from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess

def aligner(in_file); 
    """Generate an MSA of the PSI-BLAST results
    Args 
        input psi_blast file
    Returns 
        Fasta file"""

        
    out_file = "aligned.fasta"

    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=Tr

    print(clustalomega_cline)

    r = subprocess.Popen(clustalomega_cline)

