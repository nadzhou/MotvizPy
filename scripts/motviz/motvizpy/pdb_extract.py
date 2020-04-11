from Bio import SeqIO    


class MotifTable: 
    """Tabulate the possible motifs"""

    def __init__(self, seq_dict): 
        self.seq_dict = seq_dict

    def find_motif_pos(self, pos, seq):
        seq = "".join(seq)
        seq_residues = []

        for x,_ in enumerate(pos): 
            seq_residues.append(seq[x : x + 3])

        return seq_residues


    def original_file_seq_extract(self, path): 
        record = SeqIO.parse(path, "fasta")

        seq = [str(i.seq) for i in record]

        return seq  