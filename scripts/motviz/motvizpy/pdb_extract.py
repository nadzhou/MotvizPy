#!/bin/bash/env python3

from Bio import SeqIO
import urllib.parse
import urllib.request
from pathlib import Path


import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt

from stats import Analysis
from stats import seq_extract

from Bio.PDB import PDBList
from Bio.PDB import *

import re

from pypdb import describe_pdb

class MotifTable:
    """Tabulate the possible motifs"""

    def retrieve_pdb_id(self, query_ids_list):
        url = 'https://www.uniprot.org/uploadlists/'

        query_ids = " ".join(query_ids_list)

        params = {
            'from': 'P_REFSEQ_AC',
            'to': 'ACC',
            'format': 'tab',
            'query': query_ids
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
        
        return no_duplicates


    def find_motif_pos(self, pos, seq, motif_stretch=3):
        """Find the sstretches of motifs with a given threshold

        Args:
            pos [list]: 1d list of the motif positions
            seq [list]: 1d list of the sequences
            motif_stretch [int]: Position stretch to get sequence

        Returns:
            seq_residues [list]: Motifs traced from the seq from pos

        """

        seq = "".join(seq)
        seq_residues = []

        for x,_ in enumerate(pos):
            seq_residues.append(seq[x : x + motif_stretch])
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


    def pdb2fasta(self, file_path): 
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
                    

    def fasta2seq(self, file_path): 
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
                    seqs[x.id] = x.seq
        return seqs


    def find_resis_in_pdb_seqs(self, seq_residues, seq_dict): 
        found_in_seqs = []
        for k, v in seq_dict.items(): 
            for i in seq_residues: 
                if i in v: 
                    found_in_seqs.append(k)

        return found_in_seqs
    

    def describe_pdb_file(self, pdb_id, seq_residues): 
        """Retrieve information from PyPDB on protein structure

        Args: 
            pdb_id [list]: List of the PDB ids

        returns: 
            imp_data [dict]: Dictionary that contains, id, expmethod, resolution 
        
        """

        info_needed = ["structureId", "expMethod", "resolution","nr_residues", "nr_atoms", "status"]
        imp_data = []

        for id, resid in zip(pdb_id,seq_residues): 
            results = describe_pdb(id)
            temp = []
            for k, v in results.items(): 
                if k in info_needed: 
                    temp.append(v)
                    temp.append(resid)
            imp_data.append(temp)

        return imp_data


    def csv_writer(self, pdb_data, file_path): 
        """Write the dictionary data as CSV file
        
        Args: 
            pdb_data [list of lists]: List of lists of PDB info on id, exptmethod etc
            file_path [str]: File path for writing the CSV file

        """
        df = pd.DataFrame.from_records(pdb_data, \
            columns=["PDB_ID", "ExpMethod", "Rsolution", "Number of residues", "Number of atoms", "Status", "Residues"], index=False)
        df.to_csv(file_path, index=False)




def main(): 
    motif_inst = MotifTable()

    in_file = "/home/nadzhou/Desktop/omega.aln"
    seq_dict = seq_extract(in_file, "clustal")

    file_path = "/home/nadzhou/Desktop/pdbs"

    seq_dict = motif_inst.fasta2seq(file_path)
    seq_keys = list(seq_dict.keys())

    orig_file = "/home/nadzhou/Desktop/1xef.fasta"
    seq2 = motif_inst.original_file_seq_extract(orig_file)

    pos = [23, 42, 59, 123, 170, 115, 12, 13, 233]
    seq_residues = motif_inst.find_motif_pos(pos, seq2)

    found_in_seqs = motif_inst.find_resis_in_pdb_seqs(seq_residues, seq_dict)    

    pdb_data = motif_inst.describe_pdb_file(found_in_seqs, seq_resi)
    print(pdb_data)
    motif_inst.csv_writer(pdb_data, "/home/nadzhou/Desktop/pdb_results.csv")




if __name__ == '__main__': 
    main()


