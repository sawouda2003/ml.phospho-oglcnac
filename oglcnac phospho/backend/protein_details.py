import requests
import re
import pandas as pd
import time
import random
import os
import logging
import json
from csv import DictReader, DictWriter
from concurrent.futures import ThreadPoolExecutor, as_completed
import certifi
from io import StringIO

# Function to get protein details and additional data
def get_protein_details(protein_accession, max_retries=5):
    retry_delay = 1  # Initial delay in seconds
    for attempt in range(max_retries):
        try:
            protein_response = requests.post(url=f"https://api.glygen.org/protein/detail/{protein_accession}", json={"uniprot_canonical_ac": protein_accession}) #adds in the unique uniprot_ac to the url
            protein_response.raise_for_status()
            return protein_response.json()
        except requests.RequestException:
            time.sleep(retry_delay)
            retry_delay *= 2  # Double the delay for the next attempt
            retry_delay += random.uniform(0, 1)  # Add jitter

    raise Exception("Maximum retry attempts reached")

# Function to get domain details and additional data
def get_domain_details(protein_accession, max_retries=5):
    retry_delay = 1  # Initial delay in seconds
    for attempt in range(max_retries):
        try:
            protein_response = requests.get(url=f"https://www.ebi.ac.uk/proteins/api/proteins/{protein_accession}", headers={'Accept': 'application/json'}) #adds in the unique uniprot_ac to the url
            protein_response.raise_for_status()
            return protein_response.json()
        except requests.RequestException:
            print(protein_accession)
            time.sleep(retry_delay)
            retry_delay *= 2  # Double the delay for the next attempt
            retry_delay += random.uniform(0, 1)  # Add jitter

    raise Exception("Maximum retry attempts reached")

# Function to find gene name and kinase gene names
def find_gene_names_and_kinases(protein_data):
    gene_name = ""
    kinase_dict = {}
    organism = ""

    if protein_data.get('gene_names'): #if gene_names exists in the data
        for g in protein_data['gene_names']: # retrieve items in gene_names
            if (g['resource'] == 'UniProtKB') and (g['type'] == 'recommended'): 
                gene_name = g['name'] #retrieve gene name
                break
            else:
                if (g['resource'] == 'UniProtKB'):
                    gene_name = g['name']
                    break

    if protein_data.get('phosphorylation'): #if phosphorylation exists in the data
        for k in protein_data['phosphorylation']: # retrieve items in phosphorylation
            start_pos = k.get('start_pos') #get the site position
            kinase_gene_name = k.get('kinase_gene_name') #get the kinase_gene_name that acts on the given site
            if start_pos is not None and kinase_gene_name is not None:
                if start_pos in kinase_dict:
                    kinase_dict[start_pos].append(kinase_gene_name) #if the site already exists in the dictionary append the kinase_gene_name into the list
                else:
                    kinase_dict[start_pos] = [kinase_gene_name] #if the site doesnt exist create a new key value pair with site as the key and a list of kinase_gene_name as values
        for start_pos in kinase_dict:
            kinase_dict[start_pos] = ';'.join(kinase_dict[start_pos]) #finally turn the list into a string where each item in the list is joined by a ; 
    
    if protein_data.get('species'): #if species exist in the data
        for o in protein_data.get('species'): #retrieve items in species
            organism = o['name'] #retrieve the organism name for those species
    return gene_name, kinase_dict, organism

def find_domain_and_range(domain_data):
    domain_dict = {}
    
    if domain_data.get('features'): #if features exist in data
        for d in domain_data['features']: #for item in features
            if d['type'] == "DOMAIN": #if type of feature is domain
                begin = d.get('begin') #get the begin amino acid site for the domain
                end = d.get('end') #get the end amino acid site for the domain
                end = end.lstrip(">") #some sites have > in the beginning, remove that 
                desc = d.get('description') #retrieve domain description
                range = f'{begin}-{end}' #format into a range
                domain_dict[range] = desc
    return domain_dict

def find_range_and_value(data, value): #retrieving the key and value stored in domain_dict
    try:
        for key, val in data.items():
            start, end = map(int, key.split('-'))
            if start <= value <= end:
                return key, val
    except Exception as e:
        return '', ''
    return '', ''

# Old function to find GO annotations
def find_go_annotation(protein_accession, organism):
    linked_csv_files = {
        'Homo sapiens': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/human_protein_go_annotation.csv',
        'Mus musculus': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/mouse_protein_go_annotation.csv',
        'Gallus gallus': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/chicken_protein_go_annotation.csv',
        'Dictyostelium discoideum': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/dicty_protein_go_annotation.csv',
        'Drosophila melanogaster': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/fruitfly_protein_go_annotation.csv',
        'Hepacivirus C genotype 1a': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/hcv1a_protein_go_annotation.csv',
        'Hepacivirus C genotype 1b': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/hcv1b_protein_go_annotation.csv',
        'Sus scrofa': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/pig_protein_go_annotation.csv',
        'Rattus norvegicus': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/rat_protein_go_annotation.csv',
        'Severe acute respiratory syndrome coronavirus 1': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/sarscov1_protein_go_annotation.csv',
        'Severe acute respiratory syndrome coronavirus 2': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/sarscov2_protein_go_annotation.csv',
        'Saccharomyces cerevisiae S288C': 'https://data.glygen.org/ln2data/releases/data/current/reviewed/yeast_protein_go_annotation.csv'
    } #dict of csv files of the gene ontology

    logger = logging.getLogger(__name__)
    
    link = linked_csv_files.get(organism)
    if not link:
        logger.info(f"No CSV file link found for organism: {organism}")
        return '', '', ''

    try:
        response = requests.get(link, verify=True)  # Ensure SSL verification
        response.raise_for_status()  # Raise an exception for bad responses

        # Use StringIO to convert the response content to a file-like object
        csv_data = response.content.decode('utf-8')
        linked_df = pd.read_csv(StringIO(csv_data))  # Read CSV from string
    except Exception as e:
        logger.info(f"Error reading CSV file for organism: {organism} from link: {link}. Error: {e}")
        return '', '', ''
    biological_list, molecular_list, cellular_list = [], [], []
    categories = ["molecular_function", "biological_process", "cellular_component"]

    try:
        matched_df = linked_df[linked_df['uniprotkb_canonical_ac'] == protein_accession]
    except Exception as e:
        logger.info(f"Error reading CSV file for organism: {organism} from link: {link}. Error: {e}")
        return '', '', ''

    matched_df.loc[:, 'go_term_id'] = matched_df['go_term_id'].astype(str).str.replace('_', ':')
    for category in categories:  # extracting info based on categories into the corresponding list
        go_terms = matched_df.loc[matched_df['go_term_category'] == category, ['go_term_label', 'go_term_id']]
        if category == "molecular_function":
            for index, row in go_terms.iterrows():
                molecular_list.append(f"{row['go_term_label']} ({row['go_term_id']})")
        elif category == "biological_process":
            for index, row in go_terms.iterrows():
                biological_list.append(f"{row['go_term_label']} ({row['go_term_id']})")
        elif category == "cellular_component":
            for index, row in go_terms.iterrows():
                cellular_list.append(f"{row['go_term_label']} ({row['go_term_id']})")

    return '; '.join(map(str, molecular_list)), '; '.join(map(str, biological_list)), '; '.join(map(str, cellular_list))

def process_protein(row, cache, logger):
    start_time = time.time() #time recorded for logging
    uniprot = row["uniprotkb_canonical_ac"] #retrieving the protein ac
    protein_ac = re.sub(r'-\d+', '', protein_ac) #editing the protein ac to be appropriate for EBI protein API
    protein_data = get_protein_details(uniprot, max_retries=5) #retrieve protein details from API
    if protein_ac:
        try:
            domain_data = get_domain_details(protein_ac, max_retries=5) #retrieve protein details from API
            domain_dict = find_domain_and_range(domain_data)
            range, domain = find_range_and_value(domain_dict, int(row['site']))
        except Exception as e:
            logger.info(f"no domain data for {uniprot}")
            domain_dict = {}
            domain = ''
            range = ''
    gene_name, kinase_dict, organism = find_gene_names_and_kinases(protein_data) #getting the gene_name, kinase_dict, organism
    molecular_function, biological_process, cellular_component = find_go_annotation(uniprot, organism) #retrieving GO annotations
    fasta_seq = protein_data["sequence"]['sequence'] #retrieving the fasta_seq

    #making the dict to add into the cache
    target_data = {
        "gene_name": gene_name,
        "kinase_dict": kinase_dict,
        "organism": organism,
        "molecular_function": molecular_function,
        "biological_process": biological_process,
        "cellular_component": cellular_component,
        "fasta_seq": fasta_seq,
        "domain_dict" : domain_dict
    }

    #adding into the rows
    row["gene_name"] = gene_name
    row["kinase_gene_name"] = kinase_dict.get(int(row['site']))
    row['organism'] = organism
    row["molecular_function"] = molecular_function
    row["biological_process"] = biological_process
    row["cellular_component"] = cellular_component
    row["range"] = range
    row["domain"] = domain


    logger.info(f"Processed {uniprot} in {time.time() - start_time} seconds")

    with open("./output/cache.json", "w") as cache_file:
        json.dump(cache, cache_file)

    return row, target_data


# Populate additional fields and cache
def get_protein_data(use_cache, logger):
    
    #if the user uses cache, access the cache for info. else start from the cache from scratch 
    if use_cache and os.path.isfile("./output/cache.json"):
         with open("./output/cache.json", "r") as cache_file:
            cache = json.load(cache_file)
    else:
        cache = {}
    
    cache_checkpoint = 100

    columns = ["uniprotkb_canonical_ac", "protein_name", "gene_name", "site", "amino_acid", "phosphorylation", "o-glcnacylation",
           "peptide_seq_five_before", "peptide_seq_five_after", "peptide_seq_ten_before", "peptide_seq_ten_after", "kinase_gene_name", "organism",
           "taxonomy_id", "molecular_function", "biological_process", "cellular_component", "domain", "range", "status"]

    with open("./output/supersearch_results.tsv", "r") as infile, open("./output/supersearch_results.tsv" + ".tmp", "w", newline='') as outfile: # open a tmp file to be written
        reader = DictReader(infile, delimiter='\t')
        writer = DictWriter(outfile, fieldnames=columns, delimiter='\t')
        writer.writeheader()

        api_hit = 0

        # Use ThreadPoolExecutor for concurrent API requests
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = []

            for i, row in enumerate(reader): 

                uniprot = row["uniprotkb_canonical_ac"]
                if uniprot in cache: #if the uniprot has already been added to cache retrieve information from it
                    target_data = cache[uniprot]
                    row["gene_name"] = target_data["gene_name"]
                    row["kinase_gene_name"] = target_data["kinase_dict"].get(int(row['site']), '')
                    row['organism'] = target_data['organism']
                    row["molecular_function"] = target_data["molecular_function"]
                    row["biological_process"] = target_data["biological_process"]
                    row["cellular_component"] = target_data["cellular_component"]
                    range, domain = find_range_and_value(target_data['domain_dict'], int(row['site']))
                    row["range"] = range
                    row["domain"] = domain
                    writer.writerow(row)
                else:
                    # Submit the task to the executor for new proteins that have not been called from the API
                    futures.append(executor.submit(process_protein, row, cache, logger))
                
                if i % cache_checkpoint == 0:
                    results = [future.result() for future in as_completed(futures)]

                    for row, target_data in results:
                        writer.writerow(row)
                        cache[row["uniprotkb_canonical_ac"]] = target_data
                        api_hit += 1

                    # Clear the futures list after processing
                    futures.clear()

                    # Save the cache to file
                    with open("./output/cache.json", "w") as cache_file:
                        json.dump(cache, cache_file)
                    
                    logger.info("Hit checkpoint at index " + str(i))
                    logger.info("API calls " + str(api_hit))
                    logger.info("cache hits " + str(cache_checkpoint - api_hit))
                    logger.info("============================================")
                
                    api_hit = 0

            if futures:
                results = [future.result() for future in as_completed(futures)]
                for row, target_data in results:
                    writer.writerow(row)
                    cache[row["uniprotkb_canonical_ac"]] = target_data
                    api_hit += 1

            # Clear futures list after final processing
            futures.clear()

    # Replace the original file with the updated file
    os.replace("./output/supersearch_results.tsv" + ".tmp", "./output/supersearch_results.tsv")

    # Save cache
    with open("./output/cache.json", "w") as cache_file:
        json.dump(cache, cache_file)

