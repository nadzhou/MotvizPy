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
    Pdb_id = args.id_input
    pdbl = PDBList()
    ppb = PPBuilder()

    pdbl.retrieve_pdb_file(Pdb_id, file_format='pdb')
    p = Path(f'yu/pdb{Pdb_id}.ent')
    p.replace(f'{Pdb_id}.pdb')
    print("PDB file written.")

def seq_extract(args): 
    pdb_file = f"{args.id_input}.pdb"
    ppb = PPBuilder()
    p = PDBParser()

    struct = p.get_structure('x', pdb_file)
    wtih open(f"{pdb_file}.fasta", "w") as f: 
        for pp in ppb.build_peptides(struct): 
            f.write(pp.get_sequence())
    print("Sequence file written.   ")

args = parse_arguments()    
struct_seq_retrieve(args)
seq_extract(args)