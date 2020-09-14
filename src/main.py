#!/usr/bin/env python3.7

#!/usr/bin/env python

from motvizpy.seq_retrieve import parse_arguments
from motvizpy.seq_retrieve import StructSeqRetrieve
from motvizpy.psiblast import psi_blaster
from motvizpy.aligner import xml_parser
from motvizpy.msa import msa

from motvizpy.read_write_seqs import extract_seq
from motvizpy.read_write_seqs import write_seq
from motvizpy.stats import Analysis
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np

from motvizpy.identical_sequence_parser import IdenticalSequencesParser
from motvizpy.emboss import emboss_water

from Bio import AlignIO
from Bio import SeqIO

from multiprocessing import Pool

class Motviz: 
    def __init__(self, args): 
        self.pdb_id = args.id_input

        self.out_dir = Path(args.output_path)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.id_score = args.idscore

        self._motifs = None
    
    def struct_extract(self): 
        pdb_inst = StructSeqRetrieve(self.pdb_id, self.out_dir)
        pdb_inst.struct_retrieve()
        pdb_inst.replace_ent2pdb()


    def pdb_seq_extract(self): 
        pdb_seq = extract_seq(self.out_dir / f"{self.pdb_id}.pdb", "pdb-seqres")
        write_seq(pdb_seq, self.out_dir / f"{self.pdb_id}.fasta")


    def psi_blast(self): 
        psi_blaster(self.out_dir / f"{self.pdb_id}.fasta", 
                                    f"{self.out_dir}/psi.xml")

    def parse_psi_xml(self): 
        return xml_parser(self.out_dir / "psi.xml")

    def water(self): 
        emboss_water(self.out_dir / f"{self.pdb_id}.fasta",
                     self.out_dir / "seqs.fasta", 
                     self.out_dir / "water.fasta")

    def seq_trimmer(self): 
        needle_record = list(AlignIO.parse(self.out_dir / "water.fasta", "msf"))
        self.result_record = []

        for rec in needle_record[1:]: 
            reference_seq = rec[0]
            seq_parser = IdenticalSequencesParser(reference_seq, rec[1], self.id_score)

            result = seq_parser.highly_identical_seqs()
            if result:
                self.result_record.append(result)


    def clustalomega(self): 
        msa(self.out_dir / "trimmed_seqs.fasta", 
            self.out_dir / "aligned_seq.fasta")
    

    def motif_finder(self): 
        seqs = extract_seq(self.out_dir / "aligned_seq.fasta", "fasta")
        seq = [[x for x in y] for y in seqs]

        c =  Analysis(seq)                                                                                                                                                                                                                                                                                                                                                           (seq)

        self.c_ent = c.conservation_score()
        self.norm_data = c.normalize_data(self.c_ent)

        minima = c.find_local_minima(self.norm_data)
        pos_motif, pos = c.find_motif(self.norm_data, minima, threshold=-0.75)

        print(f"Positions: {pos}")
        print(f"Possible motifs: {pos_motif}")

        c.pymol_script_writer(self.out_dir / f"{self.pdb_id}_pymol.txt", pos)


    def plotter(self): 
        norm_data_len = np.arange(1, len(self.c_ent) + 1)

        sns.lineplot(norm_data_len, self.norm_data)
        plt.title(f"Conservation score per amino acid position for PDB ID {self.pdb_id}")
        plt.xlabel("Amino acid position")
        plt.ylabel("Normalized conservation score")
        plt.show() 


    def motifs(self): 
        if self._motifs == None: 
            self.struct_extract()
            self.pdb_seq_extract()
            print("PDB structure and sequence fetching: done.")

            self.psi_blast()
            print("\nParsing PSI-BLAST results...")
            psi_results = self.parse_psi_xml()
            write_seq(psi_results, self.out_dir / "seqs.fasta")

            self.water()
            self.seq_trimmer()

            try:  
                write_seq(self.result_record, 
                        self.out_dir / "trimmed_seqs.fasta")

                self.clustalomega()
                self.motif_finder()
                self.plotter()
                plt.show()

            except IOError as e: 
                print("Could not find any meaningful target sequences.\nMotif finding failed: ", e)


def main(): 
    args = parse_arguments()
    mot = Motviz(args)
    mot.motifs()


if __name__ == '__main__': 
    main()
