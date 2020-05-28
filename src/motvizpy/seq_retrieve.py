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
    
    def __init__(self, args, out_directory): 
        """Initialize the class StructSeqRetrieve
        
        Args: 
            args [argparse object]: PDB ID arg from the terminal
            out_directory [str]: File path where the files should be written
        """
        
        self.args = args
        self.out_dir = out_directory
    
    def struct_retrieve(self): 
        """Retrieve the PDB structure from the terminal argument. 
        
        Args: 
            args [argparse object]: Contains the id_input argument
            
        Returns: 
            prompt [str]: File successfully written
            
        """
        
        pdb_id = self.args.id_input
        pdbl = PDBList()
        ppb = PPBuilder()
        
        path = Path(self.out_dir)

        pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=path)
        p = Path(f"{path}/pdb{pdb_id}.ent")
        p.replace(f'{path}/{pdb_id}.pdb')
        
        return ("PDB file written.")

    def seq_extract(self): 
        """Extract the sequence from the given PDB file 
        and write it to a fasta file. 
        
        Args: 
            args [argparse object]: Take in the PDB ID
            
        Returns: 
            prompt [str]: File written to directory. 
            
        """
        
        pdb_id = f"{self.out_dir}/{self.args.id_input}.pdb"

        pdb_record = SeqIO.parse(pdb_id, "pdb-seqres")

        SeqIO.write(pdb_record, f"{self.out_dir}/{self.args.id_input}.fasta", "fasta")

        print("Fasta file written")