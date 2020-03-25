#!/usr/bin/env python3.7

from motvizpy.seq_retrieve import parse_arguments
from motvizpy.seq_retrieve import StructSeqRetrieve
from motvizpy.psiblast import psi_blaster
from motvizpy.aligner import xml_parser
from motvizpy.msa import msa
from motvizpy.stats import seq_extract
from motvizpy.stats import Analysis

pdb_inst = parse_arguments()
pdb_inst = StructSeqRetrieve(pdb_id)
pdb_inst.struct_retrieve()
pdb_inst.seq_extract()


