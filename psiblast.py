from Bio.Blast.Applications import NcbipsiblastCommandline
import subprocess



def psi_blaster(in_file):

    cmd = (f"psiblast -query {in_file} -db refseq_protein.00 -out out_psi.xml -outfmt 5")


psi_blaster()