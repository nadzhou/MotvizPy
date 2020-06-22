from aligner import xml_parser



def test_xml_parser(): 
    test_id = "ref|YP_007188579.1|"

    seqs = xml_parser("test.xml")
    assert seqs[1].id == test_id


