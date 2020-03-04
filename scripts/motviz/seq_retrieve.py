import argparse as ap 
from Bio.PDB import PDBList
from Bio.PDB import *
from pathlib import Path
from Bio import SeqIO


def parse_arguments(parser=None): 
    """Parse arguments given by the terminal for PDB ID."""
    if not parser: 
        parser = ap.ArgumentParser()
    parser.add_argument("id_input", help="Input PDB ID")
    args = parser.parse_args()

    return args

def struct_seq_retrieve(args): 
    """
    Retrieve PDB structure given the PDB Id.
    
    Args
        PDB ID
    Returns 
        PDB structure file
    """
    Pdb_id = args.id_input
    pdbl = PDBList()
    ppb = PPBuilder()

    pdbl.retrieve_pdb_file(Pdb_id, file_format='pdb', pdir=".")
    p = Path(f'pdb{Pdb_id}.ent')
    p.replace(f'{Pdb_id}.pdb')
    print("PDB file written.")

def seq_extract(args): 
    """
    Once the PDB structure is retrieved, retrieve sequence.
    Args: 
        PDB ID input
    Returns
        Sequence fasta file
    """
    pdb_file = f"{args.id_input}.pdb"
    ppb = PPBuilder()
    p = PDBParser()

    struct = p.get_structure('x', pdb_file)

    with open (f"{args.id_input}.fasta", "w") as f: 
    # OPen a fasta file and write the extracted sequence
        for pp in ppb.build_peptides(struct): 
            seq = str(pp.get_sequence())
            f.write("\n".join(seq[i : i + 80] \
                        for i in range(0, len(seq), 80)))

    print("Sequence file written.   ")

args = parse_arguments()    
struct_seq_retrieve(args)
seq_extract(args)