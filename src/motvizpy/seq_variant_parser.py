from Bio import AlignIO
import numpy as np 
import pandas as pd

from Bio import SeqIO
from pdb_extract import MotifTable
from seq_retrieve import StructSeqRetrieve


def main(): 
    seq_dir = "/home/nadzhou/SEQs/tb/needle2.fasta"
    seq_variant_parser(seq_dir)


def seq_variant_parser(seq_dir: str): 
    needle_record = list(AlignIO.parse(seq_dir, "msf"))

    for rec in needle_record: 
        consecutive_align(rec[0], rec[1])
        

def consecutive_align(seq_a: str, seq_b: str): 
    aligned_results = _align_consecutive(seq_a.seq, seq_b.seq)

    for result in aligned_results: 
        inst = MotifTable()
        map_respnose = inst.retrieve_pdb_id(seq_a.id)

        if map_respnose: 
            for ids in map_respnose: 
                retrieve_struct(ids, "/home/nadzhou/SEQs/pdbs/")
            print(seq_a.id)
            print(result)


def _align_consecutive(seq_a: str, seq_b: str): 
    consecutive_match = 0 
    potential_motif = []
    motif = ""

    for start_idx, (first_seq, second_seq) in enumerate(zip(seq_a, seq_b)): 
        hamming_dist = 0 

        match = first_seq == second_seq and first_seq != "-"

        if match: 
            consecutive_match += 1
            
        else: consecutive_match = 0 

        if consecutive_match > 3: 
            # motif += seq_a[start_idx - 2 : start_idx]     
    
            for end_idx, (nested_first_seq, nested_second_seq) in \
                            enumerate(zip(seq_a, seq_b), start=start_idx): 
                nested_match = nested_first_seq == nested_second_seq and nested_first_seq != "-"
                
                if not nested_match: 
                    hamming_dist += 1 

                if hamming_dist < 4: 
                    motif += nested_first_seq

            if motif in seq_a: 
                yield motif


def retrieve_struct(pdb_id, output_dir): 
    struct_inst = StructSeqRetrieve(pdb_id, output_dir)

    struct_inst.struct_retrieve()
    struct_inst.replace_ent2pdb()


# def find_seq_in_pdb(): 


if __name__ == "__main__":
    main()




# def alt_align(seq_a, seq_b): 
#     for index,_ in enumerate(seq_a): 
#         longest_common_seq = [x[0] == x[1] for x in zip(seq_a[index : ], seq_b[index : ]) if x[0] != '-']

#         if True not in longest_common_seq or False not in longest_common_seq: 
#             pass
#         else: 
#             print()
#             print(seq_a[index : longest_common_seq.index(True)])

