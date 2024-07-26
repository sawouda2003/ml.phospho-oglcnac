import requests
import pandas as pd
import re
from functools import reduce
import time
import random
import os
import logging
import json
from csv import DictReader, DictWriter

def get_positions_excluding_sites(fasta_seq, df, protein_accession):
    # Get the specific sites to exclude
    sites_to_exclude = df.loc[df.uniprotkb_canonical_ac == protein_accession, 'site'].astype(int).values
    
    # Find positions of 'S' excluding specific sites
    s_res = [m.start() for m in re.finditer('S', fasta_seq) if (m.start() + 1) not in sites_to_exclude]
    
    # Find positions of 'T' excluding specific sites
    t_res = [m.start() for m in re.finditer('T', fasta_seq) if (m.start() + 1) not in sites_to_exclude]
    
    return s_res, t_res

def sites_data(logger, output, df):
    # Get unique uniprotkb_canonical_ac values
    unique_ac = df['uniprotkb_canonical_ac'].unique()
    unique_ac_list = unique_ac.tolist()
    print(len(unique_ac_list))

    # Load the cache from the JSON file
    with open("./output/cache.json", "r") as cache_file:
        cache = json.load(cache_file)

    count = 1
    start_time = time.time()

    for protein_accession in cache.keys():
        if count % 1000:
            logger.info("elapse time " + str(time.time() - start_time))
            logger.info("Hit checkpoint at index " + str(count))
            logger.info("============================================")
            start_time = time.time()
            
        target_data = cache[protein_accession]
        fasta_seq = target_data['fasta_seq']
        gene_name = target_data['gene_name']
        organism = target_data['organism']
        biological_process = target_data['biological_process']
        molecular_function = target_data['molecular_function']
        cellular_component = target_data['cellular_component']

        s_res, t_res = get_positions_excluding_sites(fasta_seq, df, protein_accession)

        protein_name = df.loc[df.uniprotkb_canonical_ac == protein_accession, 'protein_name'].values[0]
        tax_id = str(df.loc[df.uniprotkb_canonical_ac == protein_accession, 'taxonomy_id'].values[0])

        for s in s_res:
            s_list = [
                protein_accession,
                protein_name,
                gene_name,  
                str(s + 1),
                "Ser",
                "",
                "",
                fasta_seq[max(0, s - 5):s],
                fasta_seq[s + 1:s + 6],
                fasta_seq[max(0, s - 10):s],
                fasta_seq[s + 1:s + 11],
                "", 
                organism,
                tax_id,
                molecular_function,  
                biological_process, 
                cellular_component,
                "",
                "",
                "unknown_glycosite"  # unknown to be glycosylated or phosphorylated
            ]
        
            s_row = "\t".join(s_list) + "\n"
            with open("./output/supersearch_results.tsv", "a") as file:
                file.write(s_row)

        for t in t_res:
            t_list = [
                protein_accession,
                protein_name,
                gene_name,
                str(t + 1),
                "Thr",
                "",
                "",
                fasta_seq[max(0, t - 5):t],
                fasta_seq[t + 1:t + 6],
                fasta_seq[max(0, t - 10):t],
                fasta_seq[t + 1:t + 11],
                '',
                organism,
                tax_id,
                molecular_function,
                biological_process,
                cellular_component,
                "",
                "",
                "unknown_glycosite"
            ]
            t_row = "\t".join(t_list) + "\n"
            with open("./output/supersearch_results.tsv", "a") as file:
                file.write(t_row)

        count += 1

    if output == "csv": #if the user chose csv convert the tsv to a csv
        df = pd.read_csv("./output/supersearch_results.tsv", delimiter="\t")

        df.to_csv("./output/supersearch_results.csv", index=False)
