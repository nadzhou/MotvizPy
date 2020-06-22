from seq_retrieve import parse_arguments
from seq_retrieve import StructSeqRetrieve
import pytest 


def test_parse_arguments(): 
    parser = parse_arguments()

    parser = parser.parse_args(['1fmw', '/home/'])

    assert parser.id_input == '1fmw'
    assert parser.output_path == '/home/'



@pytest.fixture
def struct_seq(): 
    struct_seq = StructSeqRetrieve('1fmw', '.')

    return struct_seq


def test_struct_retrieve(struct_seq): 
    struct_seq.struct_retrieve()

    with open("pdb1fmw.ent", "r") as file: 
        data = file.read()

        assert "1FMW" in data


def test_replace_ent2pdb(struct_seq ): 
    struct_seq.replace_ent2pdb()

    with open("1fmw.pdb", "r") as file: 
        data = file.read()

        assert "PDB" in data







