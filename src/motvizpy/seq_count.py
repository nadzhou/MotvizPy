from Bio import SeqIO


record = list(SeqIO.parse("/home/nadzhou/SEQs/pseudomonas/translated_b.fasta", "fasta"))


print(len(record))