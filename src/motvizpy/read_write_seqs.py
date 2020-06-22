from Bio import SeqIO
from Bio import SeqRecord

from typing import List



def write_seq(record: List, out_file: str) -> SeqRecord: 
    return SeqIO.write(record, out_file, "fasta")


def extract_seq(in_file: str, file_ext: str) -> List:        
    record = SeqIO.parse(in_file, file_ext)    

    return record   