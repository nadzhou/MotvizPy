from Bio import SeqIO
import urllib.parse
import urllib.request
from pathlib import Path


import seaborn as sns 
import matplotlib.pyplot as plt

from stats import Analysis
from stats import seq_extract

from Bio.PDB import PDBList
from Bio.PDB import *

import re

class MotifTable:
    """Tabulate the possible motifs"""

    def retrieve_pdb_id(self, query_id):
        url = 'https://www.uniprot.org/uploadlists/'


        params = {
            'from': 'ID',
            'to': 'PDB_ID',
            'format': 'tab',
            'query': query_id
            }
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        response = response.decode("utf-8")
        results = " ".join(re.split('\t|\n', response))
        results = re.findall(r'[a-zA-Z0-9]{4}', results)
        no_duplicates = set(results)
        
        for i in no_duplicates: 
            print(i, end='\t')


    def find_motif_pos(self, pos, seq):
        """Find the sstretches of motifs with a given threshold

        Args:
            pos [list]: 1d list of the motif positions
            seq [list]: 1d list of the sequences

        Returns:
            seq_residues [list]: Motifs traced from the seq from pos

        """

        seq = "".join(seq)
        seq_residues = []

        for x,_ in enumerate(pos):
            seq_residues.append(seq[x : x + 4])

        return seq_residues


    def original_file_seq_extract(self, path):
        """Extract sequence of the original file to extract mtoif

        Args:
            path [Path object]: Address of the PDB fasta file

        Returns:
            seq [list]: Seuqence of the PDB structure

        """
        record = SeqIO.parse(path, "fasta")

        seq = [str(i.seq) for i in record]

        return seq

def pdb2fasta(file_path): 
    """Retrieve the fasta file, description and
     seq from irectory containing PDB files

    Args: 
        file_path [str]: Address for folder containing PDB files
    
    """
    file_path = Path(file_path)

    for file in file_path.iterdir(): 
        file_without_ext = str(file).split(".")
        file_without_ext = file_without_ext[0]

        record = SeqIO.parse(file, "pdb-seqres")

        SeqIO.write(record, f"{file_without_ext}.fasta", "fasta")
                

def fasta2seq(file_path): 
    """Extract all the fasta sequences from a directory 
    
    Args: 
        file_path [str]: Address of the directory containing fasta files

    Returns: 
        seq [dict]: Sequences with their IDs as keys and sequences as values
    
    """
    file_path = Path(file_path)
    seqs = {}

    for file in file_path.iterdir(): 
        if "fasta" in file.name: 
            record = SeqIO.parse(file, "fasta")
            for x in record: 
                seqs[x.description] = x.seq
    return seqs



def main(): 
    in_file = "/home/nadzhou/Desktop/omega.aln"
    seq_dict = seq_extract(in_file, "clustal")

    #seq_keys = list(seq_dict.keys())
    #seq_keys = " ".join(seq_keys)


    # seq2 = pdb_inst.original_file_seq_extract(orig_file)


    #seq_residues = pdb_inst.find_motif_pos(pos,  seq2)    

    # fasta_file_path = "/home/nadzhou/Desktop/pdbs"
    # c = fasta2seq(fasta_file_path)

    # struct = list(c.keys()).split(" ")

    # for i in struct: 
    #     print(i)
    # c = "".join(str(c).strip())

    # print("Start \t Pattern \t Parent sequence location")

    # for i in seq_residues: 
    #     for m in re.finditer(str(i), c): 
    #         print(f"{m.start()} \t {m.group()} \t {c[m.start() - 3 : m.end() + 5]}")

    # plotter = {}

    # for i in seq_residues: 
    #     for k, v in seq_dict.items(): 
    #         if i in v: 
    #             plotter[k] = i

    # print("ID - motif")
    # for k, v in plotter.items(): 
    #     print(f"{k, v}", end="\t\t")

    # print("\nRefseq ID - PDB")
    seq_keys = "P04637"
    pdb_inst = MotifTable()
    pdb_inst.retrieve_pdb_id(seq_keys)



if __name__ == '__main__': 
    main()