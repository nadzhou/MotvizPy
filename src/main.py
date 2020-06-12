#!/usr/bin/env python3.7

from motvizpy.seq_retrieve import parse_arguments
from motvizpy.seq_retrieve import StructSeqRetrieve
from motvizpy.psiblast import psi_blaster
from motvizpy.aligner import xml_parser
from motvizpy.msa import msa

from motvizpy.stats import seq_extract
from motvizpy.stats import Analysis
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np


def plotter(norm_data): 
    norm_data_len = np.arange(1, len(norm_data) + 1)

    fig = sns.lineplot(norm_data_len, norm_data)
    plt.title("Conservation score per amino acid position")
    plt.xlabel("Amino acid position")
    plt.ylabel("Normalized conservation score")
    return fig


def main(): 
    path = Path(("/home/nadzhou/Desktop/6x2c"))

    pdb_id = parse_arguments()
    # pdb_inst = StructSeqRetrieve(pdb_id, path)
    # print(pdb_inst.struct_retrieve())
    # print(pdb_inst.seq_extract())

    # psi_blaster(f"{path}/{pdb_id.id_input}.fasta", f"{path}/psi.xml")

    # xml_parser(f"{path}/psi.xml", f"{path}/seqs.fasta")
    # msa(f"{path}/seqs.fasta", f"{path}/aligned_seq.fasta")

    seq = seq_extract("/home/nadzhou/Desktop/6x2c/aligned_seq.fasta", "fasta")
    seq = [[x for x in y] for y in seq]

    c = Analysis(seq)
    np_seq = c.seq2np()

    c_ent = c.conservation_score(np_seq)
    norm_data = c.normalize_data(c_ent)

    print(norm_data)

    fig = plotter(norm_data)

    plt.show()

# c.pymol_script_writer(f"{path}/mol.txt", pos)




if __name__ == '__main__': 
    main()
