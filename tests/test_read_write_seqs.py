from read_write_seqs import write_seq
from read_write_seqs import extract_seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq



def test_extract_seq(): 
    results = extract_seq("test.fasta", "fasta")


    assert "MNPIHDRTSDYHKYL" in results[0]


def test_write_seqs(): 
    data = SeqRecord(Seq("MNPIHDRTSDYHKYL"), 
                    id="1FMW", 
                    description="Test data", 
                    name="!FMW test"
                    )

    assert write_seq(data, "stk")

