import requests, sys
import re
import json



def retrieve_pdbs_from_api(uniprot_query_id): 
    """Retrieve the given UniProt ID's related structures from the InterPro API

    Args: 
        uniprot_query_id [str]: UniProt ID

    Returns: 
        record [json]: result json object ready to be parsed
     
    """

    api_src = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/"

    request_id = f"taxonomy/structure/pdb/protein/uniprot/{uniprot_query_id}"

    complete_link = api_src + request_id

    r = requests.get(complete_link, headers={ "Content-Type" : "application/json"})


    if not r.ok: 
        r.raise_for_status()
        sys.exit()

    record = r.json()

    retturn record




def parse_api_results(json_record): 
    """Parse the JSON result retrieved from the InterPro API

    Args: 
        json_record [json]: JSON object containing related information on a protein

    Returns: 
        founds_list [list]: List of PDB IDs.
        
    
    """

    metadata_list = [i for i in json_record["results"] if "metadata" in i]

    found_ids = []

    for item in metadata_list: 
        if item["structure_subset"]: 
            for i in item["structure_subset"]: 
                found_ids.append(i["accession"])


    found_ids_set = set(found_ids)
    founds_list = list(found_ids_set)

    return founds_list




def main(): 





if __name__ == '__main__': 
    main()







