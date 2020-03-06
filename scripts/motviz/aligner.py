from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def xml_parser(in_file): 
    """
    The PSI-BLAST results then need to 
    get polished so as to get the fasta 
    sequence. 

    Args: 
        Input XML file 
    Returns:    
        Psi-BLAST fasta file
    """
    seqs = []
    blast_records = SearchIO.read(in_file, "blast-xml")
    for hsp in blast_records.hsps: 
        print(hsp.hit)
        seqs.append(hsp.hit)
        
    SeqIO.write(seqs, "navid.fasta", "fasta")


xml_parser("out_psi.xml")
    
