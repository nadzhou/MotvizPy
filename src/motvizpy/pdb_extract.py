#!/bin/bash/env python3

from Bio import SeqIO
import urllib.parse
import urllib.request
from pathlib import Path


import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt


import re


class MotifTable:
    """Tabulate the possible motifs"""

    def retrieve_pdb_id(self, query_ids: str):
        url = 'https://www.uniprot.org/uploadlists/'


        params = {
            'from': 'ID',
            'to': 'PDB_ID',
            'format': 'tab',
            'query': query_ids
            }

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()

        response = response.decode("utf-8")
        results = re.findall(r'\b[\d\w]{4}\b',response)

        for result in results: 
            if result != "From": 
                yield result


if __name__ == '__main__': 
    main()


