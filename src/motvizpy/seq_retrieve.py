import argparse as ap 
from Bio.PDB import PDBList
from Bio.PDB import *
from pathlib import Path
from Bio import SeqIO


def main(): 
    parser = parse_arguments()

    struct_seq = StructSeqRetrieve(parser.id_input, parser.output_path)

    struct_seq.struct_retrieve()
    struct_seq.replace_ent2pdb()
    struct_seq.seq_extract()


def parse_arguments(parser=None): 
    """Parser will take PDB ID input from terminal to retrieve the file
    
    Args: 
        parser [argparse]: Take an argparse from terminal
        
    Returns: 
        args [argparse]: This will be used for retrieveing PDB
        
    """
    
    if not parser: 
        parser = ap.ArgumentParser()

    parser.add_argument("id_input", 
                        help="Input PDB ID")

    parser.add_argument("output_path", 
                        help="Path to output directory")

    parser.add_argument('--idscore', 
                        type=float, 
                        default=0.7, 
                        help="Input identity score for trimming of the sequences")

    return parser.parse_args()


class StructSeqRetrieve: 
    """Retrieve both PDB structure and then its sequence """
    
    def __init__(self, pdb_id, out_directory): 
        """Initialize the class StructSeqRetrieve
        
        Args: 
            args [argparse object]: PDB ID arg from the terminal
            out_directory [str]: File path where the files should be written
        """
        
        self.pdb_id = pdb_id.lower()
        self.out_dir = Path(out_directory)
        print(self.out_dir)
    
    def struct_retrieve(self): 
        """Retrieve the PDB structure from the terminal argument. 
        
        Args: 
            args [argparse object]: Contains the id_input argument
            
        Returns: 
            prompt [str]: File successfully written
            
        """
        
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(self.pdb_id, file_format='pdb', pdir=self.out_dir)


    def replace_ent2pdb(self): 
        p = self.out_dir / f"pdb{self.pdb_id}.ent"
        p.replace(self.out_dir / f'{self.pdb_id}.pdb')
        

if __name__ == '__main__':
    main()