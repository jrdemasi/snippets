#!/usr/bin/env python3

import sys
import time
from Bio import Entrez


Entrez.email = "jonathan.demasi@colorado.edu"
# We should apply for an API key so we get more queries/sec
Entrez.api_key = None

def get_complete_rsids():
    rsidlist = []
    numresults = 0
    retstart = 0
    search_string = "snp_pubmed_cited[sb]"
    search_results = Entrez.read(Entrez.esearch(db="snp", term=search_string,
                                                retmax=100000, retstart=retstart, usehistory="y"))
    print("Found a total of " +
          search_results["Count"] + " results using search string '" + search_string + "'")
    numresults = search_results["Count"]
    rsidlist = rsidlist + search_results["IdList"]
    additional_queries = int(int(numresults) / 100000)
    while additional_queries != 0:
        retstart = retstart + 100000
        search_results = Entrez.read(Entrez.esearch(db="snp", term=search_string,
                                                    retmax=100000, retstart=retstart, usehistory="y"))
        rsidlist = rsidlist + search_results["IdList"]
        additional_queries = additional_queries - 1
    return(rsidlist)

def main():
    listy = get_complete_rsids()
    return()

if __name__ == '__main__':
    main()

