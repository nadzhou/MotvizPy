 from Bio.Blast.Applications import NcbipsiblastCommandline
import subprocess



def psi_blaster(): 
    cline = NcbipsiblastCommandline(help=True)

    psi_command = cline.split("")
    psi_command_list = list(psi_command)

    r = subprocess.Popen(psi_command_list)


