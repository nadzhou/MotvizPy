#!/usr/bin/env python3.7

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
from motvizpy.emboss import emboss_needle

from Bio import AlignIO
from Bio import SeqIO


def main(): 
    args = parse_arguments()
    pdb_id = args.id_input
    out_dir = Path(args.output_path)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    struct_seq(pdb_id, out_dir)

    psi_blaster(f"{out_dir}/{args.id_input}.fasta", f"{out_dir}/psi.xml")

    psi_results = xml_parser(f"{out_dir}/psi.xml")
    write_seq(psi_results, f"{out_dir}/seqs.fasta")

    emboss_needle(f"{out_dir}/{pdb_id}.fasta",
                     out_dir/"seqs.fasta", 
                     f"{out_dir}/needle.fasta")

    seq_trimmer(f"{out_dir}/needle.fasta",f"{out_dir}/trimmed_seqs.fasta")
    msa(f"{out_dir}/trimmed_seqs.fasta", f"{out_dir}/aligned_seq.fasta")
    
    motif_finder(out_dir)


def plotter(norm_data): 
    norm_data_len = np.arange(1, len(norm_data) + 1)

    fig = sns.lineplot(norm_data_len, norm_data)
    plt.title("Conservation score per amino acid position")
    plt.xlabel("Amino acid position")
    plt.ylabel("Normalized conservation score")
    return fig


def seq_trimmer(in_file, out_file): 
    needle_record = list(AlignIO.parse(in_file, "msf"))

    result_record = []

    for rec in needle_record[1:]: 
        reference_seq = rec[0]

        seq_parser = IdenticalSequencesParser(reference_seq, rec[1])


        result = seq_parser.highly_identical_seqs()

        if result:
            result_record.append(result)

        SeqIO.write(result_record, out_file, "fasta")


def motif_finder(out_dir): 
    seqs = extract_seq(f"{out_dir}/aligned_seq.fasta", "fasta")
    seq = [[x for x in y] for y in seqs]

    c = Analysis(seq)
    np_seq = c.seq2np()

    c_ent = c.conservation_score(np_seq)
    norm_data = c.normalize_data(c_ent)

    minima = c.find_local_minima(norm_data)
    pos_motif, pos = c.find_motif(norm_data, minima, threshold=-0.75)

    fig = plotter(norm_data)

    plt.show()


def struct_seq(pdb_id, out_dir): 
    pdb_inst = StructSeqRetrieve(pdb_id, out_dir)
    pdb_inst.struct_retrieve()
    pdb_inst.replace_ent2pdb()
    pdb_seq = extract_seq(out_dir / f"{pdb_id}.pdb", "pdb-seqres")
    write_seq(pdb_seq, out_dir / f"{pdb_id}.fasta")



if __name__ == '__main__': 
    main()
