#!/usr/bin/env python3.7

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def xml_parser(in_file, out_file): 
    """Parse the XML file and then give 
    out top 100 sequences into a fasta file. 

    Args: 
        in_file [str]: Path to the input file. 
        
    Returns: 
        navid.fasta [file]: Top 100 sequences written into file
    """
    
    seqs = []
    blast_records = SearchIO.read(in_file, "blast-xml")
    for hsp in blast_records.hsps: 
        print(hsp.hit)
        seqs.append(hsp.hit)
          
    return SeqIO.write(seqs[:100], out_file, "fasta")


    
