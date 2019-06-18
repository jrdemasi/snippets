#!/usr/bin/env python3

import sys
import time
from Bio import Entrez


DEBUG = True
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

def get_pmids(interm):
    # This is obsolete now, essentially, but
    # allows a user to pass a single string
    # which can be nice.
    if isinstance(interm, str):
        interm = "rs" + interm + " AND pubmed_snp_cited[sb]"
        search_results = Entrez.read(Entrez.esearch(db="pubmed", term=interm,
                                                    retmax=100000,
                                                    usehistory="y"))
        print("Found a total of " +
              search_results["Count"] + " results using search string '" + interm + "'")
        return search_results

    elif isinstance(interm, list):
        searchstring = " OR ".join(interm)
        searchstring = "(" + searchstring + ") AND pubmed_snp_cited[sb]"
        search_results = Entrez.read(Entrez.esearch(db="pubmed",
                                                    term=searchstring,
                                                    retmax=100000,
                                                    usehistory="y"))
        print("Found a total of " +
              search_results["Count"] + " results using search string'" + searchstring + "'")
    return(search_results)

def main():
    rsids = get_complete_rsids()
    if DEBUG:
        for x in rsids:
            print(x)
    for x in rsids:
        get_pmids(x)
        time.sleep(1)
    return()

if __name__ == '__main__':
    main()

