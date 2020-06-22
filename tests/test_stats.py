import pytest
import numpy as np
from stats import Analysis



@pytest.fixture
def analyse(): 
    data = Analysis(["RSLYYD", "RSLYYD"])

    return data


def test_seq2np(analyse): 
    np_seq = analyse.seq2np()
    test_np = np.asarray(["RSLYYD", "RSLYYD"], dtype='S1')
    
    assert (np_seq == test_np).all()


def test_conservation_score(analyse): 
    ent_list = analyse.conservation_score()
    expected_values = -0.0

    assert ent_list == expected_values


def test_pymol_script_writer(analyse):
    analyse.pymol_script_writer("test_pymol.txt", [2, 4])

    with open("test_pymol.txt", "r") as file: 
        data = file.read()

    assert "2" in data