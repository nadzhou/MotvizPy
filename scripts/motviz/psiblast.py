from Bio.Blast.Applications import NcbipsiblastCommandline
import subprocess

def psi_blaster(in_file):
    """
    Via the refseq_protein database present on PC, 
    PSI-BLAST is achieved via the command line. 
    
    Args: 
        input file 
    Returns 
        PSI XML File 

    This XML file will then be parsed and only 
    sequences will be written for a round of 
    MSA. 
    """
    in_file = str(in_file)
    print("Initiating PSI-BLAST...")

    cline = NcbipsiblastCommandline(query = in_file, outfmt = 5, 
                                db = "refseq_protein.00", 
                                out="/home/nadzhou/Desktop/out_psi.xml")

    cmd = str(cline)
    cmd = cmd.split(" ")
    cmd = list(cmd)
    r = subprocess.Popen(cmd)
    r.communicate()

    print("PSI-BLAST done. File written out_psi.xml")
        
psi_blaster("./1xef.fasta")