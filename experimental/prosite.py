from Bio import SeqIO
from Bio.ExPASy import ScanProsite





aligned_record = SeqIO.parse("/home/nadzhou/SEQs/spike_uniprot.fasta", "fasta")

start_end = []

for record in aligned_record: 
    prosite_handle = ScanProsite.scan(record.seq)

    prosite_result = ScanProsite.read(prosite_handle)

    for rec in prosite_result: 
        start_end.append((rec['start'], rec['stop']))

        print(record.seq[rec['start'] : rec['stop']])
        print()

print(start_end)


