import argparse as ap 
from Bio.PDB import PDBList
from Bio.PDB import *
from pathlib import Path


def parse_arguments(parser=None): 
    if not parser: 
        parser = ap.ArgumentParser()
    parser.add_argument("id_input", help="Input PDB ID")
    args = parser.parse_args()

    return args

def struct_seq_retrieve(args): 
    struct = args.id_input
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(struct, file_format='pdb')
    p = Path(f'yu/pdb{struct}.ent')
    p.replace(f'{struct}.pdb')

args = parse_arguments()    
struct_seq_retrieve(args)