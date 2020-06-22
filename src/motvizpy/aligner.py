#!/usr/bin/env python3.7

from Bio import SearchIO
from typing import List

def xml_parser(in_file: str) -> List: 
    """Parse the XML file and then give 
    out top 100 sequences into a fasta file. 

    Args: 
        in_file [str]: Path to the input file. 
        
    Returns: 
        navid.fasta [file]: Top 100 sequences written into file
    """
    
    seqs = []
    blast_records = SearchIO.parse(in_file, "blast-xml")

    for records in blast_records: 
        for record in records: 
            for hsp in record.hsps: 

                seqs.append(hsp.hit)    
                
    return seqs
