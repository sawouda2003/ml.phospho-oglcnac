import requests
import pandas as pd
import sys
import re
from functools import reduce
import time
import random
import os
import logging
from csv import DictReader, DictWriter

SEARCH_URL = "https://{}api.glygen.org/supersearch/search"
LIST_URL = "https://{}api.glygen.org/supersearch/list"


SEARCH_PAYLOAD = {
  "concept_query_list": [
    {
      "concept": "site",
      "query": {
        "aggregator": "$or",
        "aggregated_list": [],
        "unaggregated_list": [
          {
            "path": "phosphorylation_flag",
            "order": 1,
            "operator": "$eq",
            "string_value": "true"
          },
          {
            "path": "glycosylation_flag",
            "order": 0,
            "operator": "$eq",
            "string_value": "true"
          }
        ]
      }
    }
  ]
}

# the payload is the query for the Supersearch search
def supersearch_search(search_url, list_url):
    response = requests.post(url=search_url, json=SEARCH_PAYLOAD)
    search_data = response.json()  # converting response to a json file
    list_id = search_data["results_summary"]["site"]["list_id"]  # obtaining the id
    list_res = requests.post(url=list_url, json={"id":list_id})
    list_data = list_res.json()
    total_length = list_data["pagination"]["total_length"] #retrieving the total_length of the json file
    return total_length, list_id


def fetch_with_backoff(list_url, list_id, offset, limit = 11000, max_retries=5):
    retry_delay = 1  # Initial delay in seconds
    list_payload = {
        "id": list_id, "offset": offset, "sort": "hit_score", "limit": limit, "order": "desc", "filters": []
    }
    for attempt in range(max_retries):
        try:
            response = requests.post(url=list_url, json=list_payload)
            response.raise_for_status()
            return response.json()
        except requests.RequestException:
            time.sleep(retry_delay)
            retry_delay *= 2  # Double the delay for the next attempt
            retry_delay += random.uniform(0, 1)  # Add jitter

    raise Exception("Maximum retry attempts reached")

def extract_data_from_item(item): #function to extract data from the filteres list
    def replace_yes_no(value):
        return 'Y' if value.lower() == 'yes' else 'N' if value.lower() == 'no' else value

    return {
        "uniprotkb_canonical_ac": item.get("uniprot_canonical_ac", ""),
        "protein_name": item.get("protein_name", ""),
        "gene_name": "",  # Placeholder for later population
        "site": item.get("start_pos", ""),
        "amino_acid": item.get("residue", ""),
        "phosphorylation": replace_yes_no(item.get("phosphorylation", "")),
        "o-glcnacylation": replace_yes_no(item.get("glycosylation", "")),
        "peptide_seq_five_before": item.get("up_seq", "")[-5:],
        "peptide_seq_five_after": item.get("down_seq", "")[:5] if len(item.get("down_seq", "")) > 1 else "",
        "peptide_seq_ten_before": item.get("up_seq", "")[-10:],
        "peptide_seq_ten_after": item.get("down_seq", "")[:10] if len(item.get("down_seq", "")) > 1 else "",
        "kinase_gene_name": "",  # Placeholder for later population
        "organism": "",
        "taxonomy_id": item.get("tax_id", ""),
        "molecular_function": "",  # Placeholder for later population
        "biological_process": "",  # Placeholder for later population
        "cellular_component": "",  # Placeholder for later population
        "domain": "", 
        "range": "", 
        "status": "known_glycosite"
    }

def get_all_list_data(use_cache, logger, prod_flag):
    if os.path.isfile("./output/supersearch_results.tsv") and use_cache:
        logger.info("found cache file")
        return
    
    search_url = SEARCH_URL.format("") if prod_flag else SEARCH_URL.format("beta-")
    list_url = LIST_URL.format("") if prod_flag else LIST_URL.format("beta-")

    logger.info("overwriting cache")
    limit = 11000
    offset = 1
    total_length, list_id = supersearch_search(search_url, list_url)
    columns = ["uniprotkb_canonical_ac", "protein_name", "gene_name", "site", "amino_acid", "phosphorylation", "o-glcnacylation",
           "peptide_seq_five_before", "peptide_seq_five_after", "peptide_seq_ten_before", "peptide_seq_ten_after", "kinase_gene_name", "organism",
           "taxonomy_id", "molecular_function", "biological_process", "cellular_component", "domain", "range", "status"]

    with open("./output/supersearch_results.tsv", 'w', newline='') as outfile:
        writer = DictWriter(outfile, fieldnames=columns, delimiter='\t')
        writer.writeheader()

    # Fetch data and write rows
    while offset < total_length:
        full_list_data = fetch_with_backoff(list_url, list_id, offset, limit, max_retries=5)
        site_list = [item for item in full_list_data["results"] if ((item["phosphorylation"] == "yes") or #filtering for phosphorylation and O-GlcNAcylation
                    (item["glycosylation"] == "yes" and item["phosphorylation"] == "no" and item["glycosylation_type"] == "o-linked|o-glcnacylation") or
                    (item["glycosylation"] == "no") or
                    (item["glycosylation"] == "yes" and item["phosphorylation"] == "yes" and item["glycosylation_type"] == "o-linked|o-glcnacylation"))]
        with open("./output/supersearch_results.tsv", 'a', newline='') as outfile: #writing data into the file
            writer = DictWriter(outfile, fieldnames=columns, delimiter='\t')
            for item in site_list:
                extracted_data = extract_data_from_item(item)
                writer.writerow(extracted_data)

        offset += limit
        logger.info("offset: " + str(offset))