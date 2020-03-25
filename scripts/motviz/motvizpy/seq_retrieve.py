#!/usr/bin/env python3.7

import argparse as ap 
from Bio.PDB import PDBList
from Bio.PDB import *
from pathlib import Path
from Bio import SeqIO


def parse_arguments(parser=None): 
    """Parser will take PDB ID input from terminal to retrieve the file
    
    Args: 
        parser [argparse]: Take an argparse from terminal
        
    Returns: 
        args [argparse]: This will be used for retrieveing PDB
        
    """
    
    if not parser: 
        parser = ap.ArgumentParser()
    parser.add_argument("id_input", help="Input PDB ID")
    args = parser.parse_args()

    return args

class StructSeqRetrieve: 
    """Retrieve both PDB structure and then its sequence """
    
    def __init__(self, args): 
        """Initialize the class StructSeqRetrieve
        
        Args: 
            args [argparse object]: PDB ID arg from the terminal
        """
        
        self.args = args
    
    def struct_retrieve(self): 
        """Retrieve the PDB structure from the terminal argument. 
        
        Args: 
            args [argparse object]: Contains the id_input argument
            
        Returns: 
            prompt [str]: File successfully written
            
        """
        
        Pdb_id = self.args.id_input
        pdbl = PDBList()
        ppb = PPBuilder()

        pdbl.retrieve_pdb_file(Pdb_id, file_format='pdb', pdir=".")
        p = Path(f'pdb{Pdb_id}.ent')
        p.replace(f'{Pdb_id}.pdb')
        
        return ("PDB file written.")

    def seq_extract(self): 
        """Extract the sequence from the given PDB file 
        and write it to a fasta file. 
        
        Args: 
            args [argparse object]: Take in the PDB ID
            
        Returns: 
            prompt [str]: File written to directory. 
            
        """
        
        pdb_file = f"{self.args.id_input}.pdb"
        ppb = PPBuilder()
        p = PDBParser()

        struct = p.get_structure('x', pdb_file)

        with open (f"{args.id_input}.fasta", "w") as f: 
        # OPen a fasta file and write the extracted sequence
            for pp in ppb.build_peptides(struct): 
                seq = str(pp.get_sequence())
                f.write("\n".join(seq[i : i + 80] \
                            for i in range(0, len(seq), 80)))

        return (f"Sequence file written at {args.id_input}.fasta")
