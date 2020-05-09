from Bio import AlignIO
from Bio import SeqIO


in_file = "/home/nadzhou/Desktop/divided/clustal.fasta"

alignment = []

records = list(AlignIO.parse(in_file, "fasta"))


for record in records: 
    print(record)



