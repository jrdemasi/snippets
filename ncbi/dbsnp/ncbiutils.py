#!/usr/bin/env python3

import requests

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"

# Example https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=snp_pubmed_cited[sb]&retmax=200000&retstart=1000&retmode=json
def db_query(**kwargs):
    args = []
    for key, value in kwargs.items():
        args.append(key+"="+str(value))
    qstring = "&".join(args)
    resp = requests.get(BASE_URL + qstring)
    if resp.status_code == 200:
        results = resp.json()
        return(results)
    else:
        print("You've encountered an error and we can't return your results")