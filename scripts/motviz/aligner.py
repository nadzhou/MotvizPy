from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def xml_parser(in_file): 
    
    seqs = []
    blast_records = SearchIO.read(in_file, "blast-xml")
    #with open("navid.fasta", "w") as f: 
    for hsp in blast_records.hsps: 
        print(hsp.hit)
        seqs.append(hsp.hit)
        
            #f.write(str(hsp.hit))
            #f.write("\n")

    SeqIO.write(seqs, "navid.fasta", "fasta")
    #SeqIO.write(seqs, "navid.fasta", "fasta")
    # for record in blast_records: 
    #     seqs.append(record)

    # for i in seqs: 
    #     print(dir(i))
    #print(dir(blast_records.hits[0][0]))
    
    #SeqIO.write(blast_records.hits[0][0], "navid.fasta", "fasta")

    # print(seqs[0])
    # SeqIO.write(seqs, "navid.fasta", "fasta")
    # for record in blast_records: 
    #     for i in range(50): 
    #         bir = record.hits[i]
    #         print(bir)
    #         seqs.append(bir)
    



xml_parser("out_psi.xml")
    
